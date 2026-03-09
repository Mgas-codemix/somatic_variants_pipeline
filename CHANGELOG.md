# Changelog

All notable changes to the Somatic Variants Pipeline are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

---

## [Unreleased]

### Added

#### Downstream R scripts — full implementation replacing stubs

- **`scripts/annotate_target.R`**: Real BED overlap annotation (`get_on_target` logic from Variants_Analysis doc).
  - Uses **GenomicRanges** when available; falls back to data.table non-equi join.
  - Sets `on_target`, `on_bait`, and `territory` (`"on_target"` / `"on_bait"` / `"off"`) columns.

- **`scripts/selection_candidates.R`**: Full `selectionOfCandidates()` + `checkColnames()` implementation.
  - Loads `falseCancerGenes.Rdata` and removes variants in false cancer gene list.
  - Filters by `ExonicFunc.refGene` (nonsynonymous, stopgain, stoploss, frameshift, splicing) and `Func.refGene` (exonic / splicing).
  - Adds `selected` column.

- **`scripts/annotate_functional_drug.R`**: Full drug and functional annotation.
  - Maps `ExonicFunc.refGene` → mutation category (`mc`): missense, nonsense, frameshift, splice_site, in_frame, silent, other.
  - Loads `ACC_actionable_164.Rdata` → marks actionable genes as `druggable = TRUE`.
  - Parses `DGIdb_interactions.tsv` → aggregates drug names per gene → `drug_interactions` column.
  - Sets `druggable = TRUE` when gene has DGIdb interaction or is in actionable list.

- **`scripts/qc_territories.R`**: Publication-quality ggplot2 plots.
  - Faceted by type ("all" vs "selected"), `scale_y_log10()` per Variants_Analysis doc.
  - ggsci colour palette; falls back to base R if ggplot2 unavailable.

- **`scripts/tmb.R`**: Real TMB computation.
  - Computes `exome_target_size_Mb` from target BED interval lengths (not hardcoded 50 Mb).
  - ggplot2 bar plot with `geom_hline()` at 1, 10, 100 mutations/Mb.
  - Writes `TMB.csv` with both "all" and "on_target" groups.

- **`scripts/maf_plots.R`**: Multi-facet MAF analysis.
  - Computes max MAF for four subsets: all, on_target, selected, selected_on_target.
  - ggplot2 point/line plot with `facet_wrap(~type)`.

- **`scripts/clinical_reports_xlsx.R`**: Styled Excel output.
  - Uses `openxlsx` (portable, no Java dependency) with bold header, blue fill, borders, auto-width columns.
  - Falls back to `xlsx` package if `openxlsx` not available.
  - Adds `FREQ_pct` (FREQ as percentage), `snp138_cosmic70` combined column.
  - Expanded OUTCOL with `mc`, `drug_interactions`, `territory`.

#### New scripts

- **`scripts/benchmark_performance.R`**: AUROC and Precision-Recall curve generation.
  - Computes ROC curve and AUC using FREQ and optionally CADD13_PHRED as classifier scores.
  - Computes PR curve and AUPRC.
  - Reads vcfeval `summary.txt` or hap.py `summary.csv` for TP/FP/FN counts.
  - Generates `benchmark_ROC.pdf`, `benchmark_PrecisionRecall.pdf`, `benchmark_summary.csv`.

- **`scripts/benchmark_convert_pretty_to_vcf.py`**: Convert `*_variants_pretty.csv` → VCF.
  - Supports optional truth TSV to add `TRUTH_LABEL` INFO field for score-based evaluation.

#### New Nextflow modules

- **`modules/somatic/benchmark_variants.nf`**: Opt-in benchmarking module.
  - Runs `rtg vcfeval` or `hap.py` when available.
  - Calls `benchmark_performance.R` for AUROC/PR curves.

#### nf-core structure

- **`workflows/somatic_variants.nf`**: Main workflow logic extracted for modularity.
- **`subworkflows/local/upstream.nf`**: FASTQ → pretty CSV upstream subworkflow.
- **`subworkflows/local/downstream.nf`**: pretty CSV → reports downstream subworkflow.
- **`lib/WorkflowMain.groovy`**: Pipeline-level validation utilities.
- **`lib/WorkflowSomatic.groovy`**: Samplesheet validation helpers.
- **`conf/modules.config`**: Per-process `publishDir` and `ext.args` configuration.
- **`conf/test.config`**: Minimal test profile.
- **`conf/test_full.config`**: Full integration test profile.
- **`conf/singularity.config`**: Singularity/Apptainer container profile.
- **`modules.json`**: nf-core module tracking stub.

#### Documentation

- `docs/usage.md`: Full usage guide (prerequisites, installation, samplesheet, parameters, troubleshooting).
- `docs/output.md`: Complete output file documentation with column schemas.
- `docs/downstream_processing.md`: Step-by-step downstream processing walkthrough with data-flow diagram.
- `docs/benchmarking.md`: Guide to performance evaluation (AUROC/PR, truth sets, tool installation).
- `README.md`: Added pipeline mermaid flowchart, badges, links to all doc pages.

#### Engineering

- `nextflow.config`: Added `singularity`, `test`, `test_full` profiles; benchmarking parameters.
- `environment.yml`: Added `r-ggplot2`, `r-ggsci`, `r-ggpubr`, `r-reshape2`, `r-openxlsx`, `r-proc`, `bioconductor-genomicranges`, `bioconductor-iranges`, `bioconductor-s4vectors`.
- `docker/Dockerfile`: Added `LABEL` metadata; updated R packages to include ggplot2, ggsci, ggpubr, reshape2, openxlsx, pROC; added Bioconductor GenomicRanges/IRanges.
- `.github/workflows/ci.yml`: CI workflow for R script syntax check and Nextflow lint.

### Changed

- `clinical_reports_xlsx.R`: Migrated Excel backend from `xlsx` (Java-dependent) to `openxlsx` (portable).
- `annotate_target.R`: Replaced placeholder ("mark all as off") with real BED overlap.
- `selection_candidates.R`: Replaced no-op placeholder with full filter logic.
- `annotate_functional_drug.R`: Replaced stub with full ExonicFunc mapping + DGIdb merge.
- `qc_territories.R`: Replaced base barplot with ggplot2 + log10 scale + facets.
- `tmb.R`: Replaced hardcoded 50 Mb with dynamic BED-based computation.
- `maf_plots.R`: Replaced single-subset barplot with multi-facet ggplot2 analysis.
