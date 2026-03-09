#!/usr/bin/env Rscript
# selectionOfCandidates + checkColnames. MVP stub: load ref Rdata if present, filter/standardize, write RDS.

args <- commandArgs(trailingOnly = TRUE)
input_rds <- gsub("^--input=", "", args[grepl("^--input=", args)])
false_cancer_genes   <- gsub("^--false_cancer_genes=", "", args[grepl("^--false_cancer_genes=", args)])
gene_categories_repo <- gsub("^--gene_categories_repo=", "", args[grepl("^--gene_categories_repo=", args)])
output <- gsub("^--output=", "", args[grepl("^--output=", args)])

df <- readRDS(input_rds)
if (file.exists(false_cancer_genes)) load(false_cancer_genes)
if (file.exists(gene_categories_repo)) load(gene_categories_repo)
# Placeholder: selectionOfCandidates logic; checkColnames
saveRDS(df, output)
message("Wrote ", output)
