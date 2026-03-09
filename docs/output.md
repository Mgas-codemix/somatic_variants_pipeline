# Output Documentation — Somatic Variants Pipeline

## Directory tree

```
results/
├── Variants/
│   ├── WXS/
│   │   └── <SAMPLE>_variants_pretty.csv        ← merged caller output (from_fastq)
│   ├── MuTect/
│   │   └── <SAMPLE>/
│   │       └── <SAMPLE>.vcf
│   └── VarScan/
│       └── <SAMPLE>/
│           └── <SAMPLE>_snv.vcf, <SAMPLE>_indel.vcf
├── Alignments/
│   ├── <SAMPLE>/                               ← trimmed/aligned BAMs
│   └── BQSR/<SAMPLE>/                          ← BQSR-recalibrated BAMs
├── Rdata/
│   ├── merged_variants.rds                     ← parsed + enriched table
│   ├── annotated_territory.rds                 ← after BED overlap annotation
│   ├── candidates.rds                          ← after selectionOfCandidates
│   └── somatic_mutations.rds                   ← final annotated table
├── Results/
│   └── Variant_calling/
│       ├── QCs/
│       │   ├── variants_on_territories.pdf     ← ggplot2 faceted bar plot
│       │   └── territory_summary.csv
│       ├── TMB_per_territories.pdf             ← TMB bar plot with reference lines
│       ├── max_MAF_per_territories.pdf         ← faceted MAF plot
│       └── max_MAF_summary.csv
│       └── CLINICAL_REPORTS/
│           ├── <SAMPLE>_freq0.05_DP10.drug.csv
│           └── <SAMPLE>_freq0.05_DP10.drug.xlsx
├── Tables/
│   └── TMB.csv                                 ← TMB values per sample / territory
└── Benchmarking/                               ← only when --run_benchmark
    └── <SAMPLE>/
        ├── benchmark_summary.csv
        ├── benchmark_ROC.pdf
        └── benchmark_PrecisionRecall.pdf
```

---

## File descriptions

### `Variants/WXS/<SAMPLE>_variants_pretty.csv`

Output of the upstream COLLECT_PRETTY step. Contains the merged output of MuTect and VarScan calls after ANNOVAR annotation.

**Key columns:** `Chr`, `Start`, `End`, `Ref`, `Alt`, `Func.refGene`, `Gene.refGene`, `ExonicFunc.refGene`, `AAChange.refGene`, `snp138`, `cosmic70`, `CLINSIG`, `CADD13_PHRED`, `DP`, `FREQ`, `method`, `id`, `n_methods`.

---

### `Rdata/somatic_mutations.rds`

The fully annotated R data frame. Contains all columns from the pretty CSV plus:

| Column | Description |
|--------|-------------|
| `on_target` | `TRUE` if variant overlaps the capture target BED |
| `on_bait` | `TRUE` if variant overlaps the bait BED |
| `territory` | `"on_target"`, `"on_bait"`, or `"off"` |
| `selected` | `TRUE` if variant passes `selectionOfCandidates` filters |
| `variant_type` | `"SNV"` or `"InDel"` |
| `mc` | Mutation category: `"missense"`, `"nonsense"`, `"frameshift"`, `"splice_site"`, `"in_frame"`, `"silent"`, `"other"` |
| `druggable` | `TRUE` if gene is in DGIdb or actionable gene list |
| `drug_interactions` | Semicolon-separated list of interacting drugs (from DGIdb) |

---

### `Results/Variant_calling/QCs/variants_on_territories.pdf`

ggplot2 faceted bar plot showing variant counts per territory (`on_target`, `on_bait`, `off`) and per sample. Two facets: **all** variants and **selected** variants (when the `selected` column is present). Y-axis on log10 scale per Variants_Analysis specification.

---

### `Tables/TMB.csv`

Tumor mutational burden table. Columns:

| Column | Description |
|--------|-------------|
| `sample` | Sample identifier |
| `territory` | `"all"` or `"on_target"` |
| `n_variants` | Count of variants in the group |
| `TMB` | Tumor mutational burden (mutations per megabase) |

TMB is computed as `n_variants / exome_target_size_Mb`, where `exome_target_size_Mb` is derived from the sum of intervals in the target BED file.

---

### `Results/Variant_calling/TMB_per_territories.pdf`

ggplot2 bar plot of TMB per sample (using on-target variants). Reference lines at 1, 10, and 100 mutations/Mb are drawn as dashed horizontal lines.

---

### `Results/Variant_calling/max_MAF_per_territories.pdf`

ggplot2 faceted point/line plot showing maximum allele frequency (FREQ) per sample for four variant subsets:
- `all` — all variants
- `on_target` — variants in target territory
- `selected` — variants passing `selectionOfCandidates`
- `selected_on_target` — selected and on-target

---

### `Results/Variant_calling/CLINICAL_REPORTS/<SAMPLE>_freq<F>_DP<D>.drug.xlsx`

Per-sample styled Excel workbook (produced with `openxlsx`). Contains only variants that pass:
1. `selected == TRUE`
2. `FREQ >= min_freq`
3. `DP >= min_dp`

Columns follow the OUTCOL definition in `clinical_reports_xlsx.R`:  
`Chr`, `Start`, `End`, `Ref`, `Alt`, `Gene.refGene`, `Func.refGene`, `ExonicFunc.refGene`, `AAChange.refGene`, `snp138_cosmic70`, `FREQ_pct`, `CADD13_PHRED`, `DP`, `variant_type`, `mc`, `druggable`, `drug_interactions`, `territory`, `sample`.

Bold header row, blue fill, and auto-width columns are applied.

---

### `Benchmarking/<SAMPLE>/benchmark_summary.csv`

Table summarising benchmarking results. Columns: `sample`, `score`, `AUROC`, `AUPRC`, `precision`, `recall`, `F1`.

### `Benchmarking/<SAMPLE>/benchmark_ROC.pdf`

ROC curve (True Positive Rate vs False Positive Rate) with AUC annotation. Generated by `scripts/benchmark_performance.R`.

### `Benchmarking/<SAMPLE>/benchmark_PrecisionRecall.pdf`

Precision-Recall curve with AUPRC annotation. See [docs/benchmarking.md](benchmarking.md) for interpretation guidance.
