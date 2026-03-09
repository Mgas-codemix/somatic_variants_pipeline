#!/usr/bin/env Rscript
# Write _variants_pretty.csv from multianno + variants table (HaTSPiL collect_annotated_variants).
# MVP: merge multianno with variants on key/id; select OUTCOL-like columns; write CSV.

args <- commandArgs(trailingOnly = TRUE)
sample_id    <- gsub("^--sample_id=", "", args[grepl("^--sample_id=", args)])
multianno    <- gsub("^--multianno=", "", args[grepl("^--multianno=", args)])
variants_csv <- gsub("^--variants_csv=", "", args[grepl("^--variants_csv=", args)])
output       <- gsub("^--output=", "", args[grepl("^--output=", args)])

if (!file.exists(multianno)) stop("Multianno file not found: ", multianno)
if (!file.exists(variants_csv)) stop("Variants CSV not found: ", variants_csv)

anno <- read.delim(multianno, check.names = FALSE)
vars <- read.csv(variants_csv, stringsAsFactors = FALSE)
anno$id <- paste0(anno$Chr, ":", anno$Start, "-", anno$End, "_", anno$Ref, "_", anno$Alt)
merged <- merge(anno, vars, by.x = "id", by.y = "key", all = FALSE)
outcol <- c("id", "Chr", "Start", "End", "Ref", "Alt", "Func.refGene", "Gene.refGene", "GeneDetail.refGene",
            "ExonicFunc.refGene", "AAChange.refGene", "snp138", "cosmic70", "CLINSIG", "CADD13_PHRED",
            "DP", "FREQ", "method")
outcol <- intersect(outcol, names(merged))
write.csv(merged[, outcol], output, row.names = FALSE)
message("Wrote ", output, " (", nrow(merged), " rows)")
