#!/usr/bin/env Rscript
# Subset selected; OUTCOL; format FREQ, CADD, snp138-cosmic70; write styled Excel per sample.
# See docs/Variants_Analysis_pipeline_mapping.md for full spec (Variants_Analysis doc).

# Report columns (OUTCOL from Variants_Analysis): standard order for clinical Excel/CSV
OUTCOL <- c("Chr", "Start", "End", "Ref", "Alt", "Gene.refGene", "Func.refGene", "ExonicFunc.refGene",
            "AAChange.refGene", "snp138_cosmic70", "FREQ_pct", "CADD13_PHRED", "DP",
            "variant_type", "mc", "druggable", "drug_interactions", "territory", "sample")

args <- commandArgs(trailingOnly = TRUE)
input_rds <- gsub("^--input=",    "", args[grepl("^--input=",    args)])
outdir    <- gsub("^--outdir=",   "", args[grepl("^--outdir=",   args)])
min_freq  <- as.numeric(gsub("^--min_freq=","", args[grepl("^--min_freq=", args)]))
min_dp    <- as.numeric(gsub("^--min_dp=",  "", args[grepl("^--min_dp=",   args)]))

if (is.na(min_freq)) min_freq <- 0.05
if (is.na(min_dp))   min_dp   <- 10

df <- readRDS(input_rds)

# Subset by selected flag when present (Variants_Analysis: subset(df, selected))
if ("selected" %in% names(df)) df <- df[df$selected %in% TRUE, , drop = FALSE]
if ("FREQ"     %in% names(df)) df <- df[suppressWarnings(as.numeric(df$FREQ)) >= min_freq, , drop = FALSE]
if ("DP"       %in% names(df)) df <- df[suppressWarnings(as.numeric(df$DP))   >= min_dp,   , drop = FALSE]
if (!"sample" %in% names(df)) df$sample <- "unknown"

# Combine snp138 and cosmic70 for report (doc: snp138-cosmic70 column)
if ("snp138" %in% names(df) || "cosmic70" %in% names(df)) {
  snp_part    <- if ("snp138"   %in% names(df)) as.character(df$snp138)   else ""
  cosmic_part <- if ("cosmic70" %in% names(df)) as.character(df$cosmic70) else ""
  df$snp138_cosmic70 <- paste(snp_part, cosmic_part, sep = " | ")
}

# FREQ as percentage string (display only)
if ("FREQ" %in% names(df)) {
  freq_num    <- suppressWarnings(as.numeric(df$FREQ))
  df$FREQ_pct <- ifelse(is.na(freq_num), NA_character_,
                        paste0(round(freq_num * 100, 2), "%"))
}

if ("CADD13_PHRED" %in% names(df))
  df$CADD13_PHRED <- round(suppressWarnings(as.numeric(df$CADD13_PHRED)), 2)

# Restrict to OUTCOL where present
outcol <- intersect(OUTCOL, names(df))
if (length(outcol) == 0L) outcol <- names(df)

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

for (s in unique(df$sample)) {
  sub     <- df[df$sample == s, outcol, drop = FALSE]
  out_csv <- file.path(outdir, paste0(s, "_freq", min_freq, "_DP", min_dp, ".drug.csv"))
  write.csv(sub, out_csv, row.names = FALSE)
  message("Wrote ", out_csv)

  # Styled Excel output using openxlsx (portable; no Java dependency)
  if (requireNamespace("openxlsx", quietly = TRUE)) {
    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, sheetName = s)

    # Header style: bold, light-blue background, thin border
    header_style <- openxlsx::createStyle(
      fontColour = "#000000", bgFill = "#BDD7EE",
      textDecoration = "bold", border = "TopBottomLeftRight",
      borderColour = "#4472C4", halign = "center"
    )
    body_style <- openxlsx::createStyle(
      border = "TopBottomLeftRight", borderColour = "#B8CCE4"
    )
    pct_style <- openxlsx::createStyle(
      numFmt = "0.00%", border = "TopBottomLeftRight", borderColour = "#B8CCE4"
    )

    openxlsx::writeData(wb, sheet = s, x = sub, startRow = 1, startCol = 1,
                        headerStyle = header_style, borders = "all",
                        borderStyle = "thin")
    # Apply body style to data rows
    if (nrow(sub) > 0L)
      openxlsx::addStyle(wb, sheet = s, style = body_style,
                         rows = seq(2, nrow(sub) + 1), cols = seq_len(ncol(sub)),
                         gridExpand = TRUE)

    # Auto-width columns
    openxlsx::setColWidths(wb, sheet = s, cols = seq_len(ncol(sub)), widths = "auto")

    out_xlsx <- sub("\\.csv$", ".xlsx", out_csv)
    openxlsx::saveWorkbook(wb, out_xlsx, overwrite = TRUE)
    message("Wrote ", out_xlsx)
  } else if (requireNamespace("xlsx", quietly = TRUE)) {
    # Fallback to xlsx (requires rJava)
    out_xlsx <- sub("\\.csv$", ".xlsx", out_csv)
    xlsx::write.xlsx(sub, out_xlsx, sheetName = s, row.names = FALSE)
    message("Wrote ", out_xlsx)
  }
}
