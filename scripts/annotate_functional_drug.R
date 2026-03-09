#!/usr/bin/env Rscript
# ExonicFunc refGene mapping; variant_type (SNV/InDel); merge DGIdb + actionable; druggable flag.
# MVP stub: add variant_type from Ref/Alt length, write somatic_mutations.rds.

args <- commandArgs(trailingOnly = TRUE)
input_rds   <- gsub("^--input=", "", args[grepl("^--input=", args)])
acc_path    <- gsub("^--acc_actionable=", "", args[grepl("^--acc_actionable=", args)])
dgidb_path  <- gsub("^--dgidb=", "", args[grepl("^--dgidb=", args)])
output      <- gsub("^--output=", "", args[grepl("^--output=", args)])

df <- readRDS(input_rds)
if ("Ref" %in% names(df) && "Alt" %in% names(df)) {
  df$variant_type <- ifelse(nchar(as.character(df$Ref)) == 1L & nchar(as.character(df$Alt)) == 1L, "SNV", "InDel")
} else df$variant_type <- NA
df$druggable <- NA
if (nzchar(dgidb_path) && file.exists(dgidb_path)) {
  dgidb <- read.delim(dgidb_path, check.names = FALSE)
}
saveRDS(df, output)
message("Wrote ", output)
