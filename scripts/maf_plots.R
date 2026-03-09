#!/usr/bin/env Rscript
# Max MAF overall and in target; plot max_MAF_per_territories.pdf. MVP stub.

args <- commandArgs(trailingOnly = TRUE)
input_rds <- gsub("^--input=", "", args[grepl("^--input=", args)])
outdir    <- gsub("^--outdir=", "", args[grepl("^--outdir=", args)])
pdf_out   <- gsub("^--pdf=", "", args[grepl("^--pdf=", args)])

df <- readRDS(input_rds)
if (!"FREQ" %in% names(df)) df$FREQ <- NA
if (!"sample" %in% names(df)) df$sample <- "unknown"
maf_df <- aggregate(FREQ ~ sample, data = df, FUN = function(x) max(as.numeric(x), na.rm = TRUE))
write.csv(maf_df, "max_MAF_summary.csv", row.names = FALSE)
pdf(pdf_out, width = 6, height = 4)
barplot(maf_df$FREQ, names.arg = maf_df$sample, las = 2, main = "Max MAF per sample")
dev.off()
message("Wrote ", pdf_out)
