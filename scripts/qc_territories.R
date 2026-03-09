#!/usr/bin/env Rscript
# Counts by territory (on_target, on_bait, off); bar plot variants_on_territories.pdf. MVP stub.

args <- commandArgs(trailingOnly = TRUE)
input_rds <- gsub("^--input=", "", args[grepl("^--input=", args)])
target_bed <- gsub("^--target_bed=", "", args[grepl("^--target_bed=", args)])
outdir     <- gsub("^--outdir=", "", args[grepl("^--outdir=", args)])
pdf_out    <- gsub("^--pdf=", "", args[grepl("^--pdf=", args)])

df <- readRDS(input_rds)
if ("territory" %in% names(df)) {
  counts <- as.data.frame(table(df$territory, df$sample, useNA = "ifany"))
  names(counts) <- c("territory", "sample", "count")
} else {
  counts <- data.frame(territory = "off", sample = if ("sample" %in% names(df)) unique(df$sample)[1] else "unknown", count = nrow(df))
}
write.csv(counts, "territory_summary.csv", row.names = FALSE)
pdf(pdf_out, width = 6, height = 4)
barplot(setNames(counts$count, paste(counts$territory, counts$sample)), las = 2, main = "Variants by territory")
dev.off()
message("Wrote ", pdf_out)
