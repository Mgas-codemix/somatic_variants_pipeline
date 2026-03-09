#!/usr/bin/env Rscript
# Subset selected; OUTCOL; format FREQ, CADD, snp138-cosmic70; write styled Excel per sample.
# See docs/Variants_Analysis_pipeline_mapping.md for full spec (Variants_Analysis doc).

# Report columns (OUTCOL from Variants_Analysis): standard order for clinical Excel/CSV
OUTCOL <- c("Chr", "Start", "End", "Ref", "Alt", "Gene.refGene", "Func.refGene", "ExonicFunc.refGene",
            "AAChange.refGene", "snp138_cosmic70", "FREQ", "CADD13_PHRED", "DP", "variant_type", "druggable", "sample")

args <- commandArgs(trailingOnly = TRUE)
input_rds <- gsub("^--input=", "", args[grepl("^--input=", args)])
outdir    <- gsub("^--outdir=", "", args[grepl("^--outdir=", args)])
min_freq  <- as.numeric(gsub("^--min_freq=", "", args[grepl("^--min_freq=", args)]))
min_dp    <- as.numeric(gsub("^--min_dp=", "", args[grepl("^--min_dp=", args)]))

df <- readRDS(input_rds)
# Subset by selected flag when present (Variants_Analysis: subset(df, selected))
if ("selected" %in% names(df)) df <- df[df$selected %in% TRUE, ]
if ("FREQ" %in% names(df)) df <- df[as.numeric(df$FREQ) >= min_freq, ]
if ("DP" %in% names(df))   df <- df[as.numeric(df$DP) >= min_dp, ]
if (!"sample" %in% names(df)) df$sample <- "unknown"
# Combine snp138 and cosmic70 for report (doc: snp138-cosmic70 column)
if ("snp138" %in% names(df) || "cosmic70" %in% names(df)) {
  df$snp138_cosmic70 <- paste(
    if ("snp138" %in% names(df)) as.character(df$snp138) else "",
    if ("cosmic70" %in% names(df)) as.character(df$cosmic70) else "",
    sep = " | "
  )
}
if ("CADD13_PHRED" %in% names(df))
  df$CADD13_PHRED <- round(as.numeric(df$CADD13_PHRED), 2)
# Restrict to OUTCOL where present
outcol <- intersect(OUTCOL, names(df))
if (length(outcol) == 0) outcol <- names(df)

for (s in unique(df$sample)) {
  sub <- df[df$sample == s, outcol, drop = FALSE]
  out_csv <- file.path(outdir, paste0(s, "_freq", min_freq, "_DP", min_dp, ".drug.csv"))
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  if (requireNamespace("xlsx", quietly = TRUE)) {
    out_xlsx <- sub("\\.csv$", ".xlsx", out_csv)
    xlsx::write.xlsx(sub, out_xlsx, sheetName = s, row.names = FALSE)
    message("Wrote ", out_xlsx)
  }
  write.csv(sub, out_csv, row.names = FALSE)
  message("Wrote ", out_csv)
}
