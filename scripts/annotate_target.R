#!/usr/bin/env Rscript
# Annotate variants with on_target / on_bait / off (get_on_target logic).
# MVP stub: reads RDS and BEDs, adds territory column, writes RDS. Full: use GenomicRanges/bedr overlap.

args <- commandArgs(trailingOnly = TRUE)
input_rds   <- gsub("^--input=", "", args[grepl("^--input=", args)])
target_bed  <- gsub("^--target_bed=", "", args[grepl("^--target_bed=", args)])
bait_bed    <- gsub("^--bait_bed=", "", args[grepl("^--bait_bed=", args)])
output      <- gsub("^--output=", "", args[grepl("^--output=", args)])

df <- readRDS(input_rds)
# Placeholder: if Chr/Start/End exist, could intersect with BED; else mark all as "off" for MVP
if (!"on_target" %in% names(df)) df$on_target <- NA
if (!"on_bait"   %in% names(df)) df$on_bait   <- NA
if (!"territory" %in% names(df)) df$territory <- "off"
saveRDS(df, output)
message("Wrote ", output)
