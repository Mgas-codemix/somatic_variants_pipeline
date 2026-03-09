#!/usr/bin/env Rscript
# TMB = counts / exome_target_Mb; plot with ref lines; export TMB table. MVP stub.

args <- commandArgs(trailingOnly = TRUE)
input_rds  <- gsub("^--input=", "", args[grepl("^--input=", args)])
target_bed <- gsub("^--target_bed=", "", args[grepl("^--target_bed=", args)])
outdir     <- gsub("^--outdir=", "", args[grepl("^--outdir=", args)])
tmb_csv    <- gsub("^--tmb_csv=", "", args[grepl("^--tmb_csv=", args)])
ref_lines  <- gsub("^--ref_lines=", "", args[grepl("^--ref_lines=", args)])

df <- readRDS(input_rds)
if (!"sample" %in% names(df)) df$sample <- "unknown"
exome_mb <- 50
tmb_df <- aggregate(list(n = rep(1, nrow(df))), by = list(sample = df$sample), length)
tmb_df$TMB <- tmb_df$n / exome_mb
write.csv(tmb_df, tmb_csv, row.names = FALSE)
pdf("TMB_per_territories.pdf", width = 6, height = 4)
barplot(tmb_df$TMB, names.arg = tmb_df$sample, las = 2, main = "TMB (per exome Mb)"); abline(h = c(1, 10, 100), lty = 2)
dev.off()
message("Wrote ", tmb_csv)
