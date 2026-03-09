#!/usr/bin/env Rscript
# benchmark_performance.R
# Reads TP/FP/FN counts from a benchmarking output (vcfeval / hap.py summary),
# sweeps FREQ thresholds to compute Precision-Recall curve,
# computes AUROC using variant quality scores vs truth labels,
# and generates publication-quality plots:
#   benchmark_ROC.pdf               — ROC curve with AUC annotation
#   benchmark_PrecisionRecall.pdf   — PR curve with AUPRC annotation
#   benchmark_summary.csv           — Table per sample/caller

args <- commandArgs(trailingOnly = TRUE)
bench_dir  <- gsub("^--bench_dir=",  "", args[grepl("^--bench_dir=",  args)])
pretty_csv <- gsub("^--pretty_csv=", "", args[grepl("^--pretty_csv=", args)])
truth_vcf  <- gsub("^--truth_vcf=",  "", args[grepl("^--truth_vcf=",  args)])
outdir     <- gsub("^--outdir=",     "", args[grepl("^--outdir=",     args)])
sample_id  <- gsub("^--sample=",     "", args[grepl("^--sample=",     args)])
if (!nzchar(sample_id)) sample_id <- "unknown"

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# Read benchmark counts (vcfeval produces summary.txt; hap.py produces summary.csv)
# ---------------------------------------------------------------------------
read_vcfeval_summary <- function(path) {
  if (!file.exists(path)) return(NULL)
  lines <- readLines(path)
  # vcfeval summary has header line and data rows with TP/FP/FN fields
  # Try to parse as whitespace-delimited table
  tryCatch({
    df <- read.table(text = paste(lines, collapse = "\n"), header = TRUE,
                     fill = TRUE, comment.char = "#", stringsAsFactors = FALSE)
    df
  }, error = function(e) NULL)
}

read_happy_summary <- function(path) {
  if (!file.exists(path)) return(NULL)
  tryCatch(read.csv(path, stringsAsFactors = FALSE), error = function(e) NULL)
}

# Try to locate benchmark results
bench_summary <- NULL
if (nzchar(bench_dir) && dir.exists(bench_dir)) {
  # vcfeval output
  vsummary <- file.path(bench_dir, "summary.txt")
  hsummary <- file.path(bench_dir, "summary.csv")
  if (file.exists(vsummary)) bench_summary <- read_vcfeval_summary(vsummary)
  if (is.null(bench_summary) && file.exists(hsummary)) bench_summary <- read_happy_summary(hsummary)
}

# ---------------------------------------------------------------------------
# Load pretty CSV for score-based ROC/PR computation
# ---------------------------------------------------------------------------
variants_df <- NULL
if (nzchar(pretty_csv) && file.exists(pretty_csv)) {
  library(data.table)
  variants_df <- as.data.frame(fread(pretty_csv, data.table = FALSE))
  message(sprintf("Loaded %d variants from %s", nrow(variants_df), pretty_csv))
}

# ---------------------------------------------------------------------------
# Helper: compute ROC and PR curves from scores + labels
# Returns list(fpr, tpr, auc, precision, recall, auprc, thresholds)
# ---------------------------------------------------------------------------
compute_curves <- function(scores, labels) {
  # labels: 1 = TP (true positive in truth set), 0 = FP
  ord       <- order(scores, decreasing = TRUE)
  scores    <- scores[ord]
  labels    <- labels[ord]
  n_pos     <- sum(labels == 1L)
  n_neg     <- sum(labels == 0L)
  if (n_pos == 0L || n_neg == 0L) return(NULL)

  tp <- cumsum(labels == 1L)
  fp <- cumsum(labels == 0L)
  tpr  <- tp / n_pos
  fpr  <- fp / n_neg
  prec <- tp / (tp + fp)
  rec  <- tpr

  # AUC (trapezoidal) — sort by fpr to ensure correct integration direction
  ord_fpr <- order(fpr)
  fpr_s   <- fpr[ord_fpr]
  tpr_s   <- tpr[ord_fpr]
  auc     <- sum(diff(fpr_s) * (tpr_s[-1] + tpr_s[-length(tpr_s)]) / 2)

  # AUPRC — sort by recall (independent variable for PR)
  ord_rec <- order(rec)
  rec_s   <- rec[ord_rec]
  prec_s  <- prec[ord_rec]
  auprc   <- sum(diff(rec_s) * (prec_s[-1] + prec_s[-length(prec_s)]) / 2)

  list(fpr = fpr_s, tpr = tpr_s, auc = abs(auc),
       precision = prec_s, recall = rec_s, auprc = abs(auprc),
       thresholds = scores[ord_fpr])
}

# ---------------------------------------------------------------------------
# Precision-Recall sweep across FREQ thresholds using bench_summary TP/FP/FN
# ---------------------------------------------------------------------------
summary_rows <- list()
curves       <- NULL

if (!is.null(variants_df) &&
    "FREQ" %in% names(variants_df) &&
    "truth_label" %in% names(variants_df)) {
  # truth_label must be present (0/1); added by benchmark_convert_pretty_to_vcf.py
  freq_vals <- suppressWarnings(as.numeric(variants_df$FREQ))
  labels    <- as.integer(variants_df$truth_label)
  valid     <- !is.na(freq_vals) & !is.na(labels) & labels %in% c(0L, 1L)
  if (sum(valid) > 5L) {
    curves <- compute_curves(freq_vals[valid], labels[valid])
    message(sprintf("ROC AUC=%.3f   PR AUPRC=%.3f", curves$auc, curves$auprc))
    summary_rows[["freq_roc"]] <- data.frame(
      sample  = sample_id, score  = "FREQ",
      AUROC   = round(curves$auc,   3),
      AUPRC   = round(curves$auprc, 3),
      stringsAsFactors = FALSE
    )
  }

  # Optional: CADD score as classifier
  if ("CADD13_PHRED" %in% names(variants_df)) {
    cadd_vals <- suppressWarnings(as.numeric(variants_df$CADD13_PHRED))
    valid2    <- !is.na(cadd_vals) & !is.na(labels) & labels %in% c(0L, 1L)
    if (sum(valid2) > 5L) {
      cadd_curves <- compute_curves(cadd_vals[valid2], labels[valid2])
      summary_rows[["cadd_roc"]] <- data.frame(
        sample  = sample_id, score  = "CADD13_PHRED",
        AUROC   = round(cadd_curves$auc,   3),
        AUPRC   = round(cadd_curves$auprc, 3),
        stringsAsFactors = FALSE
      )
    }
  }
}

# Fallback: compute precision/recall from bench_summary counts
if (!is.null(bench_summary)) {
  for (col_tp in c("TP", "True.Positives.Truth", "TP_base")) {
    if (col_tp %in% names(bench_summary)) {
      tp  <- suppressWarnings(as.numeric(bench_summary[[col_tp]][1]))
      fp  <- suppressWarnings(as.numeric(bench_summary[[intersect(c("FP","False.Positives","FP_call"), names(bench_summary))[1]]][1]))
      fn  <- suppressWarnings(as.numeric(bench_summary[[intersect(c("FN","False.Negatives","FN_base"), names(bench_summary))[1]]][1]))
      if (!is.na(tp) && !is.na(fp) && !is.na(fn) && (tp + fp) > 0 && (tp + fn) > 0) {
        prec <- tp / (tp + fp)
        rec  <- tp / (tp + fn)
        f1   <- 2 * prec * rec / (prec + rec)
        summary_rows[["bench"]] <- data.frame(
          sample    = sample_id, score = "vcfeval",
          precision = round(prec, 3), recall = round(rec, 3),
          F1        = round(f1,   3), AUROC  = NA, AUPRC = NA,
          stringsAsFactors = FALSE
        )
      }
      break
    }
  }
}

# Write summary CSV
summary_df <- do.call(rbind, lapply(summary_rows, function(x) {
  for (col in c("precision","recall","F1","AUROC","AUPRC"))
    if (!col %in% names(x)) x[[col]] <- NA
  x
}))
if (!is.null(summary_df) && nrow(summary_df) > 0L) {
  write.csv(summary_df, file.path(outdir, "benchmark_summary.csv"), row.names = FALSE)
  message("Wrote benchmark_summary.csv")
}

# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------
use_gg <- requireNamespace("ggplot2", quietly = TRUE)

if (!is.null(curves)) {
  # ---- ROC ----
  roc_df <- data.frame(fpr = curves$fpr, tpr = curves$tpr)

  if (use_gg) {
    library(ggplot2)
    p_roc <- ggplot(roc_df, aes(x = fpr, y = tpr)) +
      geom_line(colour = "#2166AC", linewidth = 1) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
      annotate("text", x = 0.7, y = 0.1,
               label = sprintf("AUC = %.3f", curves$auc), size = 4) +
      labs(title = sprintf("ROC Curve — %s", sample_id),
           x = "False Positive Rate", y = "True Positive Rate") +
      theme_bw(base_size = 12)
    pdf(file.path(outdir, "benchmark_ROC.pdf"), width = 6, height = 6)
    print(p_roc)
    dev.off()
    message("Wrote benchmark_ROC.pdf")
  } else {
    pdf(file.path(outdir, "benchmark_ROC.pdf"), width = 6, height = 6)
    plot(curves$fpr, curves$tpr, type = "l", col = "#2166AC", lwd = 2,
         xlab = "FPR", ylab = "TPR", main = sprintf("ROC (AUC=%.3f)", curves$auc))
    abline(0, 1, lty = 2, col = "grey60")
    dev.off()
    message("Wrote benchmark_ROC.pdf")
  }

  # ---- Precision-Recall ----
  pr_df <- data.frame(recall = curves$recall, precision = curves$precision)

  if (use_gg) {
    p_pr <- ggplot(pr_df, aes(x = recall, y = precision)) +
      geom_line(colour = "#D6604D", linewidth = 1) +
      annotate("text", x = 0.2, y = 0.1,
               label = sprintf("AUPRC = %.3f", curves$auprc), size = 4) +
      labs(title = sprintf("Precision-Recall Curve — %s", sample_id),
           x = "Recall", y = "Precision") +
      coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
      theme_bw(base_size = 12)
    pdf(file.path(outdir, "benchmark_PrecisionRecall.pdf"), width = 6, height = 6)
    print(p_pr)
    dev.off()
    message("Wrote benchmark_PrecisionRecall.pdf")
  } else {
    pdf(file.path(outdir, "benchmark_PrecisionRecall.pdf"), width = 6, height = 6)
    plot(curves$recall, curves$precision, type = "l", col = "#D6604D", lwd = 2,
         xlab = "Recall", ylab = "Precision",
         main = sprintf("PR curve (AUPRC=%.3f)", curves$auprc),
         xlim = c(0, 1), ylim = c(0, 1))
    dev.off()
    message("Wrote benchmark_PrecisionRecall.pdf")
  }
} else {
  message("No score/label data available for curve generation. ",
          "Provide a pretty CSV with a 'truth_label' column.")
}

message("benchmark_performance.R complete.")
