#!/usr/bin/env Rscript
# TMB = variant counts / exome_target_size_Mb.
# Computes exome_target_size_Mb from target BED (sum of interval lengths).
# Writes TMB.csv and TMB_per_territories.pdf.

args <- commandArgs(trailingOnly = TRUE)
input_rds  <- gsub("^--input=",     "", args[grepl("^--input=",     args)])
target_bed <- gsub("^--target_bed=","", args[grepl("^--target_bed=",args)])
outdir     <- gsub("^--outdir=",    "", args[grepl("^--outdir=",    args)])
tmb_csv    <- gsub("^--tmb_csv=",   "", args[grepl("^--tmb_csv=",   args)])
ref_lines  <- gsub("^--ref_lines=", "", args[grepl("^--ref_lines=", args)])

library(data.table)

# ---- Compute exome target size from BED ----
exome_mb <- 50  # fallback: ~50 Mb is a typical human whole-exome capture size; used only when target BED is missing/unreadable
if (nzchar(target_bed) && file.exists(target_bed)) {
  bed <- tryCatch(
    fread(target_bed, header = FALSE, sep = "\t", select = 1:3,
          col.names = c("chr", "start", "end"), data.table = TRUE),
    error = function(e) NULL
  )
  if (!is.null(bed) && nrow(bed) > 0L) {
    bed[, len := as.integer(end) - as.integer(start)]
    exome_mb <- sum(bed$len, na.rm = TRUE) / 1e6
    message(sprintf("Computed exome target size: %.3f Mb from %s", exome_mb, target_bed))
  }
} else {
  message(sprintf("target_bed not found – using default exome size of %g Mb", exome_mb))
}

# ---- Load data ----
df <- readRDS(input_rds)
if (!"sample"    %in% names(df)) df$sample    <- "unknown"
if (!"territory" %in% names(df)) df$territory <- "off"

# ---- Overall TMB (all on-target variants per sample) ----
on_target_df <- df[df$territory == "on_target", , drop = FALSE]
total_df     <- df

count_tmb <- function(sub, label) {
  ct <- as.data.frame(table(sample = sub$sample), stringsAsFactors = FALSE)
  names(ct)[names(ct) == "Freq"] <- "n_variants"
  ct$territory  <- label
  ct$TMB        <- ct$n_variants / exome_mb
  ct
}

tmb_df <- rbind(
  count_tmb(total_df,     "all"),
  count_tmb(on_target_df, "on_target")
)

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
write.csv(tmb_df, tmb_csv, row.names = FALSE)
message("Wrote ", tmb_csv)

# ---- Reference lines ----
ref_vals <- tryCatch(as.numeric(strsplit(ref_lines, "[,;]+")[[1]]),
                     error = function(e) c(1, 10, 100))
if (length(ref_vals) == 0L || all(is.na(ref_vals))) ref_vals <- c(1, 10, 100)

# ---- Plot ----
use_gg <- requireNamespace("ggplot2", quietly = TRUE)

if (use_gg) {
  library(ggplot2)
  tmb_plot <- tmb_df[tmb_df$territory == "on_target", , drop = FALSE]
  if (nrow(tmb_plot) == 0L) tmb_plot <- tmb_df[tmb_df$territory == "all", , drop = FALSE]

  p <- ggplot(tmb_plot, aes(x = sample, y = TMB, fill = sample)) +
    geom_col() +
    geom_hline(yintercept = ref_vals, linetype = "dashed", colour = "grey40") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(title = sprintf("Tumor Mutational Burden (exome: %.2f Mb)", exome_mb),
         x     = "Sample",
         y     = "TMB (mutations / Mb)",
         fill  = "Sample") +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  if (requireNamespace("ggsci", quietly = TRUE))
    p <- p + ggsci::scale_fill_npg()

  pdf_path <- file.path(outdir, "TMB_per_territories.pdf")
  pdf(pdf_path, width = 7, height = 5)
  print(p)
  dev.off()
  message("Wrote ", pdf_path)
} else {
  on_tgt <- tmb_df[tmb_df$territory == "on_target", , drop = FALSE]
  if (nrow(on_tgt) == 0L) on_tgt <- tmb_df
  pdf_path <- file.path(outdir, "TMB_per_territories.pdf")
  pdf(pdf_path, width = 6, height = 4)
  barplot(on_tgt$TMB, names.arg = on_tgt$sample, las = 2,
          main = sprintf("TMB (exome: %.2f Mb)", exome_mb),
          ylab = "TMB (mutations/Mb)")
  abline(h = ref_vals, lty = 2, col = "grey40")
  dev.off()
  message("Wrote ", pdf_path)
}
