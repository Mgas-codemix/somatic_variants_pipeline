#!/usr/bin/env Rscript
# Parse and merge *_variants_pretty.csv; compute n_methods, DP.A, add sample column.
# MVP stub: reads CSVs, adds sample_id, writes RDS. Full logic: n_methods from method column, DP.A = FREQ * DP.

library(data.table)
args <- commandArgs(trailingOnly = TRUE)
csv_paths <- strsplit(gsub("^--csv_paths=", "", args[grepl("^--csv_paths=", args)]), ",")[[1]]
sample_ids <- strsplit(gsub("^--sample_ids=", "", args[grepl("^--sample_ids=", args)]), ",")[[1]]
output    <- gsub("^--output=", "", args[grepl("^--output=", args)])

if (length(csv_paths) != length(sample_ids)) stop("csv_paths and sample_ids length mismatch")
out_list <- list()
for (i in seq_along(csv_paths)) {
  p <- csv_paths[i]
  if (!file.exists(p)) { warning("File not found: ", p); next }
  d <- tryCatch(fread(p), error = function(e) read.csv(p, sep = "\t", check.names = FALSE))
  d$sample <- sample_ids[i]
  if ("method" %in% names(d)) d$n_methods <- lengths(strsplit(as.character(d$method), ":"))
  if (all(c("FREQ", "DP") %in% names(d))) d$DP.A <- as.numeric(d$FREQ) * as.numeric(d$DP)
  out_list[[i]] <- d
}
if (length(out_list) == 0) stop("No valid CSV files")
df <- rbindlist(out_list, fill = TRUE)
saveRDS(df, output)
message("Wrote ", output, " (", nrow(df), " rows)")
