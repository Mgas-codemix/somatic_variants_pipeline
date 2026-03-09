#!/usr/bin/env Rscript
# Counts by territory (on_target, on_bait, off); ggplot2 bar plot (faceted by type:
# "all" variants and "selected" variants). Saves variants_on_territories.pdf.

args <- commandArgs(trailingOnly = TRUE)
input_rds  <- gsub("^--input=",     "", args[grepl("^--input=",     args)])
target_bed <- gsub("^--target_bed=","", args[grepl("^--target_bed=",args)])
outdir     <- gsub("^--outdir=",    "", args[grepl("^--outdir=",    args)])
pdf_out    <- gsub("^--pdf=",       "", args[grepl("^--pdf=",       args)])

library(data.table)

df <- readRDS(input_rds)
if (!"sample"    %in% names(df)) df$sample    <- "unknown"
if (!"territory" %in% names(df)) df$territory <- "off"

make_counts <- function(sub, type_label) {
  ct <- as.data.frame(table(territory = sub$territory,
                             sample    = sub$sample,
                             useNA     = "ifany"),
                      stringsAsFactors = FALSE)
  ct$type <- type_label
  ct
}

counts_all <- make_counts(df, "all")
counts <- counts_all

if ("selected" %in% names(df)) {
  sel_df      <- df[!is.na(df$selected) & df$selected %in% TRUE, , drop = FALSE]
  counts_sel  <- make_counts(sel_df, "selected")
  counts      <- rbind(counts_all, counts_sel)
}

names(counts)[names(counts) == "Freq"] <- "count"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
write.csv(counts, file.path(outdir, "territory_summary.csv"), row.names = FALSE)

# ---- ggplot2 plot (falls back to base if unavailable) ----
use_gg <- requireNamespace("ggplot2", quietly = TRUE)

if (use_gg) {
  library(ggplot2)
  TERR_LEVELS <- c("on_target", "on_bait", "off")
  counts$territory <- factor(counts$territory, levels = TERR_LEVELS)
  counts$type      <- factor(counts$type,      levels = c("all", "selected"))

  p <- ggplot(counts, aes(x = territory, y = count + 1, fill = sample)) +
    geom_col(position = "dodge") +
    facet_wrap(~type) +
    scale_y_log10() +
    labs(title  = "Variants by territory",
         x      = "Territory",
         y      = "Count (log10 scale)",
         fill   = "Sample") +
    theme_bw(base_size = 11)

  # Optional ggsci palette
  if (requireNamespace("ggsci", quietly = TRUE))
    p <- p + ggsci::scale_fill_npg()

  pdf(pdf_out, width = 8, height = 5)
  print(p)
  dev.off()
} else {
  # Base fallback
  pdf(pdf_out, width = 6, height = 4)
  mat <- tapply(counts_all$count,
                list(counts_all$territory, counts_all$sample),
                FUN = sum, default = 0)
  barplot(mat, beside = TRUE, las = 2, log = "y",
          main = "Variants by territory", ylab = "Count (log10)")
  dev.off()
}

message("Wrote ", pdf_out)
