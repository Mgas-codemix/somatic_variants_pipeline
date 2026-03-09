# Test Data for Somatic Variants Pipeline

This directory contains synthetic test data for end-to-end pipeline testing, covering both
`from_fastq` (upstream) and `from_pretty_csv` (downstream) entry points.

---

## Directory Contents

### FASTQ / Alignment Test Data (`from_fastq`)

| File | Description |
|------|-------------|
| `mock_tumor_R1.fastq.gz` | Synthetic paired-end FASTQ, read 1 (5 reads, 50 bp) |
| `mock_tumor_R2.fastq.gz` | Synthetic paired-end FASTQ, read 2 (5 reads, 50 bp) |
| `samplesheet_fastq.csv` | Samplesheet with one sample (`mock_tumor`, `from_fastq`) |
| `reference/ref.fa` | Minimal synthetic reference FASTA |
| `params_from_fastq_mock.config` | Nextflow config for `from_fastq` mock run |

### Downstream Variant Analysis Test Data (`from_pretty_csv`)

| File | Description |
|------|-------------|
| `sample1_variants_pretty.csv` | ~40 synthetic variants — lung cancer profile (TP53, KRAS, EGFR, BRAF + passengers) |
| `sample2_variants_pretty.csv` | ~25 synthetic variants — melanoma profile (BRAF V600E, ALK, EGFR druggable targets) |
| `sample3_variants_pretty.csv` | ~18 synthetic variants — low-TMB colorectal profile |
| `samplesheet_pretty.csv` | Samplesheet for `from_pretty_csv` with all 3 samples |
| `target.bed` | Target capture intervals overlapping test variant positions |
| `bait.bed` | Bait capture intervals (padded target) |
| `DGIdb_interactions.tsv` | Drug-gene interaction table (EGFR→Erlotinib, BRAF→Vemurafenib, etc.) |
| `falseCancerGenes.Rdata` | Character vector of passenger/false-driver genes (TTN, MUC16, …) |
| `gene_categories_repo.Rdata` | data.frame mapping genes to TSG/oncogene categories |
| `ACC_actionable_164.Rdata` | Character vector of FDA-actionable gene names |
| `generate_test_rdata.R` | R script to (re)generate the three `.Rdata` files above |
| `params_from_pretty_csv_test.config` | Nextflow config for `from_pretty_csv` test run |

### Benchmarking / Performance Evaluation Data

| File | Description |
|------|-------------|
| `truth_variants.csv` | Truth set with TP/FP/FN labels for AUROC/PR benchmarking |
| `confident_regions.bed` | High-confidence genomic regions for truth-set evaluation |

---

## Generating Reference Rdata Files

The `.Rdata` files must be generated once before running downstream tests.
From the **pipeline root** (parent of `test_data/`):

```bash
Rscript test_data/generate_test_rdata.R
```

This creates:
- `test_data/falseCancerGenes.Rdata`
- `test_data/gene_categories_repo.Rdata`
- `test_data/ACC_actionable_164.Rdata`

---

## Running the `from_fastq` Test

```bash
nextflow run main.nf -entry from_fastq \
  --samplesheet test_data/samplesheet_fastq.csv \
  --target_bed  test_data/target.bed \
  --bait_bed    test_data/bait.bed \
  --outdir      results_test_fastq \
  -profile local
```

> The `from_fastq` path requires reference data and external tools. The run will stop with a
> clear error at the first missing dependency (e.g. `--ref_fasta`, `--annovar_basedir`).

---

## Running the `from_pretty_csv` Test

Use the convenience script from the pipeline root:

```bash
./run_test_pretty.sh
```

Or directly:

```bash
# Generate Rdata files first (only needed once)
Rscript test_data/generate_test_rdata.R

nextflow run main.nf -entry from_pretty_csv \
  -c test_data/params_from_pretty_csv_test.config \
  -profile local
```

### Expected Outputs

After a successful `from_pretty_csv` run, the `results_test_pretty/` directory should contain:

| Output | Produced by |
|--------|-------------|
| `*/somatic_mutations.rds` | `ANNOTATE_FUNCTIONAL_DRUG` |
| `*/TMB_per_territories.pdf` + `*_TMB.csv` | `TMB` |
| `*/QC_territories.pdf` | `QC_TERRITORIES` |
| `*/MAF_plots.pdf` | `MAF_PLOTS` |
| `*/clinical_report_*.xlsx` | `CLINICAL_REPORTS_XLSX` |

### What Each Script Should Produce with the Test Data

| Script | Expected behaviour |
|--------|--------------------|
| `annotate_target.R` | Mix of `on_target`, `on_bait`, and `off` territories (variants overlap BEDs) |
| `selection_candidates.R` | TTN/MUC16 variants filtered out; exonic/nonsynonymous variants marked `selected = TRUE` |
| `annotate_functional_drug.R` | EGFR, BRAF, ALK, KRAS variants flagged as `druggable = TRUE` |
| `tmb.R` | Per-sample TMB computed from BED size (not hardcoded 50 Mb) |
| `qc_territories.R` | Multi-facet PDF with "all" and "selected" panels |
| `maf_plots.R` | MAF distribution plots (all / selected / on-target) |
| `clinical_reports_xlsx.R` | Per-sample Excel/CSV reports |

---

## Benchmarking the Pipeline

The `truth_variants.csv` file provides a labelled truth set for performance evaluation:

```
variant_id,Chr,Start,End,Ref,Alt,Gene,truth_label,in_pipeline_output,notes
chr17:7577120-7577120_C_T,...,TP53,TP,TRUE,...
```

Labels:
- `TP` — true positive (variant is real and pipeline should detect it)
- `FP` — false positive (variant present in output but should be filtered)
- `FN` — false negative (variant is real but not detected by pipeline)

The `confident_regions.bed` defines high-confidence genomic intervals matching `target.bed`,
ensuring truth calls outside these regions are excluded from evaluation.

After running the benchmarking module (`scripts/benchmark_performance.R`), AUROC and
Precision-Recall curves are generated per sample, per caller, and per variant type.

---

## Variant Design Principles

The synthetic pretty CSV files are designed to exercise every downstream script:

1. **Positional consistency** — variant positions overlap BED intervals so `annotate_target.R`
   produces a genuine mix of `on_target` / `on_bait` / `off` territories.
2. **Gene diversity** — driver genes (TP53, KRAS), actionable genes (EGFR, BRAF, ALK), and
   passenger/false-driver genes (TTN, MUC16).
3. **Method diversity** — `Mutect1.17` only, `VarScan2` only, and `Mutect1.17:VarScan2` consensus.
4. **Frequency spectrum** — sub-clonal (0.01-0.10), clonal (0.20-0.50), and near-homozygous (0.88-0.90).
5. **Depth range** — low-coverage (DP 5-15), standard (50-200), and deep (300-500).
6. **Variant types** — SNVs, frameshift deletions, non-frameshift deletions, stopgains, splicing.
7. **Functional diversity** — exonic/nonsynonymous, synonymous, intronic, UTR3, UTR5, splicing, intergenic.
8. **CADD scores** — low (0-5 for passengers) to high (30-40 for known hotspots).
