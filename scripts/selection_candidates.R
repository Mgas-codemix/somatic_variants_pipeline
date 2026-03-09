#!/usr/bin/env Rscript
# selectionOfCandidates + checkColnames.
# Loads falseCancerGenes.Rdata and gene_categories_repo.Rdata,
# validates required columns, filters variants, and adds a 'selected' column.

args <- commandArgs(trailingOnly = TRUE)
input_rds            <- gsub("^--input=",              "", args[grepl("^--input=",              args)])
false_cancer_genes   <- gsub("^--false_cancer_genes=", "", args[grepl("^--false_cancer_genes=", args)])
gene_categories_repo <- gsub("^--gene_categories_repo=", "", args[grepl("^--gene_categories_repo=", args)])
output               <- gsub("^--output=",             "", args[grepl("^--output=",             args)])

# ---------------------------------------------------------------------------
# checkColnames: warn (not abort) about missing columns so partial data still works
# ---------------------------------------------------------------------------
REQUIRED_COLS <- c("Chr", "Start", "End", "Ref", "Alt",
                   "Gene.refGene", "Func.refGene", "ExonicFunc.refGene",
                   "FREQ", "DP")

checkColnames <- function(df) {
  missing <- setdiff(REQUIRED_COLS, names(df))
  if (length(missing) > 0L)
    warning("Missing expected columns (will proceed): ", paste(missing, collapse = ", "))
  invisible(missing)
}

# ---------------------------------------------------------------------------
# selectionOfCandidates: apply biological filters per Variants_Analysis doc
# ---------------------------------------------------------------------------
selectionOfCandidates <- function(df, false_genes_vec = character(0)) {

  n_start <- nrow(df)

  # 1. Remove variants in false cancer gene list
  if (length(false_genes_vec) > 0L && "Gene.refGene" %in% names(df)) {
    keep <- !(df$Gene.refGene %in% false_genes_vec)
    df   <- df[keep, , drop = FALSE]
    message(sprintf("  Removed %d variants in falseCancerGenes list",
                    n_start - nrow(df)))
  }

  # 2. Filter by ExonicFunc.refGene (keep coding-relevant consequences)
  KEEP_EXONIC <- c(
    "nonsynonymous SNV",
    "stopgain",
    "stoploss",
    "frameshift deletion",
    "frameshift insertion",
    "nonframeshift deletion",
    "nonframeshift insertion",
    "splicing"
  )
  if ("ExonicFunc.refGene" %in% names(df)) {
    exonic_ok <- df$ExonicFunc.refGene %in% KEEP_EXONIC
    # Variants with no ExonicFunc (intronic splicing etc.) are handled by Func filter below
    exonic_na <- is.na(df$ExonicFunc.refGene) | df$ExonicFunc.refGene == "" |
                 df$ExonicFunc.refGene == "."
    df$selected <- exonic_ok | exonic_na   # provisional; refined by Func filter
  } else {
    df$selected <- TRUE
  }

  # 3. Filter by Func.refGene (keep exonic / splicing regions)
  KEEP_FUNC <- c("exonic", "splicing", "exonic;splicing")
  if ("Func.refGene" %in% names(df)) {
    func_ok     <- df$Func.refGene %in% KEEP_FUNC
    df$selected <- df$selected & func_ok
  }

  n_selected <- sum(df$selected, na.rm = TRUE)
  message(sprintf("  selectionOfCandidates: %d / %d variants selected",
                  n_selected, nrow(df)))
  df
}

# ---------------------------------------------------------------------------
# Load reference data
# ---------------------------------------------------------------------------
false_genes_vec <- character(0)

if (nzchar(false_cancer_genes) && file.exists(false_cancer_genes)) {
  env_fcg <- new.env(parent = emptyenv())
  load(false_cancer_genes, envir = env_fcg)
  # Look for a character/data.frame object named 'falseCancerGenes' or similar
  fcg_obj <- ls(env_fcg)
  for (nm in fcg_obj) {
    obj <- get(nm, envir = env_fcg)
    if (is.character(obj))          { false_genes_vec <- obj; break }
    if (is.data.frame(obj) && ncol(obj) >= 1L) {
      false_genes_vec <- as.character(obj[[1]]); break
    }
  }
  message(sprintf("Loaded %d false cancer genes from %s",
                  length(false_genes_vec), false_cancer_genes))
} else {
  message("false_cancer_genes not provided or not found – skipping false-gene filter")
}

if (nzchar(gene_categories_repo) && file.exists(gene_categories_repo)) {
  env_gcr <- new.env(parent = emptyenv())
  load(gene_categories_repo, envir = env_gcr)
  message("Loaded gene_categories_repo from ", gene_categories_repo)
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
df <- readRDS(input_rds)
checkColnames(df)
df <- selectionOfCandidates(df, false_genes_vec = false_genes_vec)

saveRDS(df, output)
message("Wrote ", output)
