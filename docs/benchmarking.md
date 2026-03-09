# Benchmarking and Performance Evaluation

This guide explains how to evaluate the performance of the somatic variants pipeline against a truth set, interpret AUROC and Precision-Recall curves, and customise the benchmarking workflow.

---

## Overview

Performance benchmarking is **opt-in** and activated by providing a truth VCF:

```bash
nextflow run main.nf \
  -entry from_pretty_csv \
  --samplesheet samplesheet.csv \
  --target_bed  /path/to/target.bed \
  --bait_bed    /path/to/bait.bed \
  --truth_vcf   /path/to/truth.vcf.gz \
  --confident_bed /path/to/confident.bed \
  --run_benchmark \
  --outdir results
```

When `--run_benchmark` is set, the `BENCHMARK_VARIANTS` module runs after variant annotation and generates:
- `benchmark_summary.csv` — precision, recall, F1, AUROC, AUPRC per sample
- `benchmark_ROC.pdf` — ROC curve with AUC annotation
- `benchmark_PrecisionRecall.pdf` — PR curve with AUPRC annotation

---

## Obtaining truth sets

### GIAB (Genome in a Bottle) — recommended

For NA12878 (HG001) hg38:
```bash
# Truth VCF
wget ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/latest/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
# Confident regions BED
wget ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/latest/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.bed
```

For other GIAB samples (HG002–HG007), see: https://www.nist.gov/programs-projects/genome-bottle

### Synthetic truth sets

For development/testing, generate a synthetic truth set from a subset of chromosome 22:
```bash
# Subset your calls to chr22
bcftools view -r chr22 calls.vcf.gz -O z -o calls_chr22.vcf.gz
# Use a known set of chr22 variants as truth (e.g. from dbSNP common variants)
```

---

## Variant comparison tools

The `BENCHMARK_VARIANTS` module tries the following tools in order:

1. **rtg vcfeval** (RTG Tools) — most accurate; requires an SDF-format reference
2. **hap.py** (Illumina/GA4GH) — widely used; requires a FASTA reference

If neither is available, the module falls back to `scripts/benchmark_performance.R` for score-based evaluation using the `truth_label` column added by `benchmark_convert_pretty_to_vcf.py`.

### Installing rtg vcfeval

```bash
# Conda
conda install -c bioconda rtg-tools
# Or download from https://www.realtimegenomics.com/products/rtg-tools

# Build SDF reference
rtg format -o ref.sdf /path/to/genome.fa
```

### Installing hap.py

```bash
pip install hap.py
# Or: conda install -c bioconda hap.py
```

---

## CSV to VCF conversion

Before running `vcfeval` or `hap.py`, convert the pretty CSV to VCF:

```bash
python scripts/benchmark_convert_pretty_to_vcf.py \
    --input  SAMPLE_variants_pretty.csv \
    --output SAMPLE_calls.vcf \
    --truth  truth_variants.tsv \     # optional: adds TRUTH_LABEL INFO field
    --sample SAMPLE_ID
```

When `--truth` is provided, each variant gets a `TRUTH_LABEL=1` (found in truth) or `TRUTH_LABEL=0` (not found) INFO field. This enables the score-based AUROC computation in `benchmark_performance.R`.

---

## Score-based evaluation (AUROC / PR curves)

When a pretty CSV has a `truth_label` column (0/1), `benchmark_performance.R` computes:

### ROC Curve (Receiver Operating Characteristic)

- **X-axis**: False Positive Rate (FPR = FP / (FP + TN))
- **Y-axis**: True Positive Rate / Recall (TPR = TP / (TP + FN))
- **AUC** (Area Under Curve): ranges 0.5 (random) to 1.0 (perfect)
- Classifier scores used: `FREQ` (allele frequency) and optionally `CADD13_PHRED`

### Precision-Recall Curve

- **X-axis**: Recall = TP / (TP + FN)
- **Y-axis**: Precision = TP / (TP + FP)
- **AUPRC** (Area Under PR Curve): ranges 0.0 to 1.0; more informative than AUC for imbalanced datasets (true variants are rare)

### Expected performance ranges

| Metric | Typical range | Notes |
|--------|---------------|-------|
| SNV AUROC | 0.85–0.98 | High with FREQ or CADD as score |
| InDel AUROC | 0.70–0.92 | InDels are harder to call |
| SNV AUPRC | 0.70–0.95 | Depends heavily on VAF threshold |
| Precision@Recall=0.9 | 0.80–0.95 | With min_freq=0.05 |

---

## Stratified benchmarking

To stratify by variant type, territory, or caller, post-process the `benchmark_summary.csv` or run `benchmark_performance.R` on subsets:

```r
# Stratify by variant_type
for (vt in c("SNV", "InDel")) {
  sub <- df[df$variant_type == vt, ]
  # run compute_curves(sub$FREQ, sub$truth_label)
}
```

---

## Outputs

| File | Description |
|------|-------------|
| `benchmark_summary.csv` | Per-sample/caller: precision, recall, F1, AUROC, AUPRC |
| `benchmark_ROC.pdf` | ROC curve(s) |
| `benchmark_PrecisionRecall.pdf` | Precision-Recall curve(s) |

All outputs are in `results/Benchmarking/<SAMPLE>/`.

---

## Troubleshooting

**Q: `vcfeval` fails with `ReferenceDataException`**  
A: The SDF reference does not match the VCF contig names. Ensure the same reference was used for variant calling and SDF creation.

**Q: Both `vcfeval` and `hap.py` are missing**  
A: The benchmark module falls back to score-based evaluation. Install one of the tools for variant-level comparison.

**Q: `truth_label` column is missing from pretty CSV**  
A: Run `benchmark_convert_pretty_to_vcf.py` with `--truth truth_variants.tsv` to add truth labels. Without this column, AUROC/PR curves cannot be generated.
