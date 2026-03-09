#!/usr/bin/env Rscript
# Max MAF analysis: computes max allele frequency for four variant subsets
# (all / selected / on-target / selected+on-target) and saves a faceted ggplot2 PDF.

args <- commandArgs(trailingOnly = TRUE)
input_rds <- gsub("^--input=",  "", args[grepl("^--input=",  args)])
outdir    <- gsub("^--outdir=", "", args[grepl("^--outdir=", args)])
pdf_out   <- gsub("^--pdf=",    "", args[grepl("^--pdf=",    args)])

df <- readRDS(input_rds)
if (!"FREQ"      %in% names(df)) df$FREQ      <- NA
if (!"sample"    %in% names(df)) df$sample    <- "unknown"
if (!"territory" %in% names(df)) df$territory <- "off"

df$FREQ_num <- suppressWarnings(as.numeric(df$FREQ))

compute_max_maf <- function(sub, type_label) {
  if (nrow(sub) == 0L)
    return(data.frame(sample = character(0), max_MAF = numeric(0), type = character(0)))
  agg <- aggregate(FREQ_num ~ sample, data = sub, FUN = function(x) max(x, na.rm = TRUE))
  names(agg)[names(agg) == "FREQ_num"] <- "max_MAF"
  agg$type <- type_label
  agg
}

subsets <- list(
  "all"      = df,
  "on_target"= df[!is.na(df$territory) & df$territory == "on_target", , drop = FALSE]
)

if ("selected" %in% names(df)) {
  sel_df <- df[!is.na(df$selected) & df$selected %in% TRUE, , drop = FALSE]
  subsets[["selected"]]            <- sel_df
  subsets[["selected_on_target"]]  <- sel_df[!is.na(sel_df$territory) &
                                               sel_df$territory == "on_target", , drop = FALSE]
}

maf_list <- mapply(compute_max_maf, subsets, names(subsets), SIMPLIFY = FALSE)
maf_df   <- do.call(rbind, maf_list)
row.names(maf_df) <- NULL

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
write.csv(maf_df, file.path(outdir, "max_MAF_summary.csv"), row.names = FALSE)

use_gg <- requireNamespace("ggplot2", quietly = TRUE)

if (use_gg) {
  library(ggplot2)
  TYPE_LEVELS <- c("all", "on_target", "selected", "selected_on_target")
  maf_df$type <- factor(maf_df$type, levels = intersect(TYPE_LEVELS, unique(maf_df$type)))

  p <- ggplot(maf_df, aes(x = sample, y = max_MAF, colour = sample, group = sample)) +
    geom_point(size = 3) +
    geom_line(aes(group = 1)) +
    facet_wrap(~type, scales = "free_y") +
    labs(title  = "Max MAF per sample and variant subset",
         x      = "Sample",
         y      = "Max allele frequency",
         colour = "Sample") +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  if (requireNamespace("ggsci", quietly = TRUE))
    p <- p + ggsci::scale_colour_npg()

  pdf(pdf_out, width = 9, height = 6)
  print(p)
  dev.off()
} else {
  pdf(pdf_out, width = 6, height = 4)
  all_df <- maf_df[maf_df$type == "all", , drop = FALSE]
  barplot(all_df$max_MAF, names.arg = all_df$sample, las = 2,
          main = "Max MAF per sample (all variants)")
  dev.off()
}

message("Wrote ", pdf_out)
