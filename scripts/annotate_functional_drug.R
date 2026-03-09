#!/usr/bin/env Rscript
# ExonicFunc.refGene mapping; variant_type (SNV/InDel); merge DGIdb + actionable; druggable flag.
# Saves somatic_mutations.rds.

args <- commandArgs(trailingOnly = TRUE)
input_rds  <- gsub("^--input=",        "", args[grepl("^--input=",        args)])
acc_path   <- gsub("^--acc_actionable=","", args[grepl("^--acc_actionable=",args)])
dgidb_path <- gsub("^--dgidb=",        "", args[grepl("^--dgidb=",        args)])
output     <- gsub("^--output=",       "", args[grepl("^--output=",       args)])

# ---------------------------------------------------------------------------
# ExonicFunc → mutation category (mc) mapping (Variants_Analysis doc)
# ---------------------------------------------------------------------------
exonic_to_mc <- function(ef) {
  ef <- as.character(ef)
  dplyr_map <- c(
    "nonsynonymous SNV"       = "missense",
    "stopgain"                = "nonsense",
    "stoploss"                = "nonsense",
    "frameshift deletion"     = "frameshift",
    "frameshift insertion"    = "frameshift",
    "nonframeshift deletion"  = "in_frame",
    "nonframeshift insertion" = "in_frame",
    "splicing"                = "splice_site",
    "synonymous SNV"          = "silent"
  )
  mc <- dplyr_map[ef]
  mc[is.na(mc)] <- "other"
  unname(mc)
}

# ---------------------------------------------------------------------------
# variant_type: SNV vs InDel
# ---------------------------------------------------------------------------
get_variant_type <- function(ref, alt) {
  ref <- as.character(ref)
  alt <- as.character(alt)
  ifelse(nchar(ref) == 1L & nchar(alt) == 1L, "SNV", "InDel")
}

# ---------------------------------------------------------------------------
# Load reference data
# ---------------------------------------------------------------------------
acc_genes <- character(0)
if (nzchar(acc_path) && file.exists(acc_path)) {
  env_acc <- new.env(parent = emptyenv())
  load(acc_path, envir = env_acc)
  for (nm in ls(env_acc)) {
    obj <- get(nm, envir = env_acc)
    if (is.character(obj)) { acc_genes <- obj; break }
    if (is.data.frame(obj) && ncol(obj) >= 1L) {
      acc_genes <- as.character(obj[[1]]); break
    }
  }
  message(sprintf("Loaded %d actionable genes from %s", length(acc_genes), acc_path))
}

dgidb_genes <- character(0)
dgidb_df    <- NULL
if (nzchar(dgidb_path) && file.exists(dgidb_path)) {
  dgidb_df <- tryCatch(
    read.delim(dgidb_path, header = TRUE, sep = "\t",
               stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) { warning("Could not read DGIdb file: ", e$message); NULL }
  )
  if (!is.null(dgidb_df)) {
    # DGIdb standard column names; tolerate minor variations
    gene_col <- intersect(c("gene_name", "gene", "Gene", "GENE"), names(dgidb_df))[1]
    drug_col <- intersect(c("drug_name", "drug", "Drug", "DRUG", "interaction_drug_name"),
                          names(dgidb_df))[1]
    if (!is.na(gene_col)) {
      dgidb_genes <- unique(as.character(dgidb_df[[gene_col]]))
      message(sprintf("Loaded DGIdb: %d unique genes with drug interactions", length(dgidb_genes)))

      # Aggregate drug names per gene (comma-separated, up to 5 per gene)
      if (!is.na(drug_col)) {
        agg <- tapply(as.character(dgidb_df[[drug_col]]),
                      as.character(dgidb_df[[gene_col]]),
                      function(x) paste(unique(x[!is.na(x) & x != ""])[seq_len(min(5L, length(x)))],
                                        collapse = "; "))
        dgidb_df <- data.frame(gene = names(agg), drug_interactions = as.character(agg),
                               stringsAsFactors = FALSE)
      } else {
        dgidb_df <- data.frame(gene = dgidb_genes, drug_interactions = NA_character_,
                               stringsAsFactors = FALSE)
      }
    }
  }
}

# ---------------------------------------------------------------------------
# Main annotation
# ---------------------------------------------------------------------------
df <- readRDS(input_rds)

# variant_type
if (all(c("Ref", "Alt") %in% names(df))) {
  df$variant_type <- get_variant_type(df$Ref, df$Alt)
} else {
  df$variant_type <- NA_character_
}

# mutation category (mc)
if ("ExonicFunc.refGene" %in% names(df)) {
  df$mc <- exonic_to_mc(df$ExonicFunc.refGene)
} else {
  df$mc <- "other"
}

# druggable flag and drug_interactions column
df$druggable        <- FALSE
df$drug_interactions <- NA_character_

if ("Gene.refGene" %in% names(df)) {
  genes <- as.character(df$Gene.refGene)

  # Mark actionable
  df$druggable <- df$druggable | (genes %in% acc_genes)

  # Merge DGIdb interactions
  if (!is.null(dgidb_df) && nrow(dgidb_df) > 0L) {
    idx <- match(genes, dgidb_df$gene)
    has_interaction <- !is.na(idx)
    df$druggable[has_interaction]         <- TRUE
    df$drug_interactions[has_interaction] <- dgidb_df$drug_interactions[idx[has_interaction]]
  }

  n_drug <- sum(df$druggable, na.rm = TRUE)
  message(sprintf("Druggable variants: %d / %d", n_drug, nrow(df)))
}

saveRDS(df, output)
message("Wrote ", output)
