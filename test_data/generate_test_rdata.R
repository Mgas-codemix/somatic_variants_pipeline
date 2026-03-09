#!/usr/bin/env Rscript
# Generate synthetic reference Rdata files for testing
# Run from the pipeline root directory:
#   Rscript test_data/generate_test_rdata.R

# False cancer genes (large genes with many passenger mutations due to their size,
# not because they are true oncogenic drivers)
falseCancerGenes <- c(
  "TTN", "MUC16", "OBSCN", "AHNAK2", "SYNE1", "SYNE2",
  "DNAH5", "PCLO", "RYR2", "RYR3", "USH2A", "DNAH11",
  "HMCN1", "FLG", "ZBTB20"
)
save(falseCancerGenes, file = "test_data/falseCancerGenes.Rdata")
message("Saved falseCancerGenes.Rdata (", length(falseCancerGenes), " genes)")

# Gene categories: tumour suppressor genes (TSG) and oncogenes
gene_categories_repo <- data.frame(
  gene = c(
    "TP53", "KRAS", "EGFR", "BRAF", "PIK3CA", "ALK", "MET",
    "ROS1", "NTRK1", "NTRK2", "NTRK3", "RET", "BRCA1", "BRCA2",
    "APC", "PTEN", "RB1", "NF1", "VHL", "WT1", "SMAD4", "CDH1",
    "NRAS", "HRAS", "IDH1", "IDH2", "FGFR1", "FGFR2", "FGFR3",
    "ERBB2"
  ),
  category = c(
    "TSG", "oncogene", "oncogene", "oncogene", "oncogene",
    "oncogene", "oncogene", "oncogene", "oncogene", "oncogene",
    "oncogene", "oncogene", "TSG", "TSG", "TSG", "TSG", "TSG",
    "TSG", "TSG", "TSG", "TSG", "TSG", "oncogene", "oncogene",
    "oncogene", "oncogene", "oncogene", "oncogene", "oncogene",
    "oncogene"
  ),
  stringsAsFactors = FALSE
)
save(gene_categories_repo, file = "test_data/gene_categories_repo.Rdata")
message("Saved gene_categories_repo.Rdata (", nrow(gene_categories_repo), " genes)")

# Actionable genes with FDA-approved companion diagnostics or targeted therapies
ACC_actionable_164 <- c(
  "EGFR", "BRAF", "ALK", "ROS1", "KRAS", "NTRK1", "NTRK2", "NTRK3",
  "MET", "RET", "ERBB2", "PIK3CA", "BRCA1", "BRCA2", "MSH2", "MLH1",
  "PMS2", "MSH6", "FGFR2", "FGFR3", "IDH1", "IDH2", "FLT3", "KIT",
  "PDGFRA", "ABL1", "NPM1", "CEBPA", "TP53", "ATM", "CHEK2", "PALB2",
  "RAD51C", "RAD51D", "NRAS", "HRAS", "MAP2K1", "MAP2K2", "NF1",
  "TSC1", "TSC2", "PTEN", "AKT1", "MTOR"
)
save(ACC_actionable_164, file = "test_data/ACC_actionable_164.Rdata")
message("Saved ACC_actionable_164.Rdata (", length(ACC_actionable_164), " genes)")

message("All test reference Rdata files generated successfully in test_data/")
