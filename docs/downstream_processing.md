# Downstream Processing Walkthrough

This document describes the **downstream** section of the somatic variants pipeline — from a collection of `*_variants_pretty.csv` files through to QC plots, TMB, MAF analysis, and clinical Excel reports.

---

## Data flow overview

```
*_variants_pretty.csv  (one per sample)
         │
         ▼
  ┌─────────────────────┐
  │  PARSE_PRETTY_CSV   │  parse_and_enrich.R
  │  merged_variants.rds│  n_methods, DP.A, sample columns added
  └────────┬────────────┘
           │
           ▼
  ┌─────────────────────────┐
  │  ANNOTATE_TARGET_BAIT   │  annotate_target.R
  │  annotated_territory.rds│  on_target, on_bait, territory columns
  └────────┬────────────────┘
           │
           ▼
  ┌─────────────────────────┐
  │  SELECT_CANDIDATES      │  selection_candidates.R
  │  candidates.rds         │  selected column (TRUE/FALSE)
  └────────┬────────────────┘
           │
           ▼
  ┌──────────────────────────────┐
  │  ANNOTATE_FUNCTIONAL_DRUG    │  annotate_functional_drug.R
  │  somatic_mutations.rds       │  mc, variant_type, druggable, drug_interactions
  └────────┬─────────────────────┘
           │
           ├──────────────────────┐
           │                      │
           ▼                      ▼
  ┌─────────────────┐    ┌─────────────────┐
  │ QC_TERRITORIES  │    │      TMB        │
  │ variants_on_    │    │ TMB.csv         │
  │ territories.pdf │    │ TMB_per_        │
  └─────────────────┘    │ territories.pdf │
                         └─────────────────┘
           │
           ▼
  ┌─────────────────────────────┐
  │  MAF_ANALYSIS               │
  │  max_MAF_per_territories.pdf│
  └──────────┬──────────────────┘
             │
             ▼
  ┌──────────────────────────────────────┐
  │  CLINICAL_REPORTS                    │
  │  <SAMPLE>_freq*_DP*.drug.xlsx / .csv │
  └──────────────────────────────────────┘
```

---

## Step-by-step description

### 1. PARSE_PRETTY_CSV → `merged_variants.rds`

**Script:** `scripts/parse_and_enrich.R`

- Reads all `*_variants_pretty.csv` files listed in the samplesheet
- Adds:
  - `n_methods`: number of callers that detected the variant (split `method` column by `:`)
  - `DP.A`: allele-specific depth = `FREQ × DP`
  - `sample`: sample identifier
- Combines all sample data frames into a single table with `rbindlist(fill=TRUE)`

**Output:** `Rdata/merged_variants.rds`

---

### 2. ANNOTATE_TARGET_BAIT → `annotated_territory.rds`

**Script:** `scripts/annotate_target.R`

Implements the `get_on_target()` logic from Variants_Analysis.doc:

- Reads target BED (`--target_bed`) and bait BED (`--bait_bed`)
- For each variant (Chr, Start, End):
  - Checks overlap with target intervals → `on_target = TRUE/FALSE`
  - Checks overlap with bait intervals → `on_bait = TRUE/FALSE`
- Derives `territory`:
  - `"on_target"` if `on_target == TRUE`
  - `"on_bait"` if `on_bait == TRUE` (and not on_target)
  - `"off"` otherwise
- Uses **GenomicRanges** when available (most accurate); falls back to **data.table non-equi join** (fast, no Bioconductor required)

**Output:** `Rdata/annotated_territory.rds`

---

### 3. SELECT_CANDIDATES → `candidates.rds`

**Script:** `scripts/selection_candidates.R`

Implements `selectionOfCandidates()` and `checkColnames()` from Variants_Analysis.doc:

**`checkColnames()`**: warns if any of these required columns are missing:
`Chr`, `Start`, `End`, `Ref`, `Alt`, `Gene.refGene`, `Func.refGene`, `ExonicFunc.refGene`, `FREQ`, `DP`

**`selectionOfCandidates()`** — three-step filter:

1. **False cancer gene filter**: removes variants where `Gene.refGene` is in the `falseCancerGenes.Rdata` list
2. **ExonicFunc filter**: keeps variants with ExonicFunc.refGene in:
   - `nonsynonymous SNV`, `stopgain`, `stoploss`, `frameshift deletion`, `frameshift insertion`, `nonframeshift deletion`, `nonframeshift insertion`, `splicing`
3. **Func.refGene filter**: keeps variants where Func.refGene is `exonic`, `splicing`, or `exonic;splicing`

Adds `selected = TRUE` for passing variants, `FALSE` otherwise.

**Customisation:**
- Replace `falseCancerGenes.Rdata` with a custom gene list (`--false_cancer_genes`)
- Adjust the `KEEP_EXONIC` / `KEEP_FUNC` vectors in `selection_candidates.R`

**Output:** `Rdata/candidates.rds`

---

### 4. ANNOTATE_FUNCTIONAL_DRUG → `somatic_mutations.rds`

**Script:** `scripts/annotate_functional_drug.R`

Adds three new annotation tracks:

#### variant_type
Derived from Ref/Alt length:
- `"SNV"` if `nchar(Ref) == 1 && nchar(Alt) == 1`
- `"InDel"` otherwise

#### mc (mutation category)
Maps `ExonicFunc.refGene` → standard category:

| ExonicFunc.refGene | mc |
|--------------------|----|
| nonsynonymous SNV | missense |
| stopgain / stoploss | nonsense |
| frameshift deletion/insertion | frameshift |
| nonframeshift deletion/insertion | in_frame |
| splicing | splice_site |
| synonymous SNV | silent |
| other | other |

#### druggable + drug_interactions
- Loads `ACC_actionable_164.Rdata` → marks genes as `druggable = TRUE`
- Parses `DGIdb_interactions.tsv` → aggregates drug names per gene → merges on `Gene.refGene`
- Sets `drug_interactions` to a semicolon-separated list of up to 5 drugs per gene

**Output:** `Rdata/somatic_mutations.rds`

---

### 5. QC_TERRITORIES → `variants_on_territories.pdf`

**Script:** `scripts/qc_territories.R`

- Counts variants by territory (`on_target`, `on_bait`, `off`) and sample
- Two facets: **"all"** variants and **"selected"** variants (when `selected` column present)
- Y-axis on `log10` scale (per Variants_Analysis doc)
- Uses **ggplot2** + **ggsci** colour palette when available; falls back to base R

**Output:** `Results/Variant_calling/QCs/variants_on_territories.pdf`, `territory_summary.csv`

---

### 6. TMB → `TMB.csv` and `TMB_per_territories.pdf`

**Script:** `scripts/tmb.R`

- Computes `exome_target_size_Mb` = sum of BED interval lengths / 1e6 (from `--target_bed`)
- Computes TMB per sample for two groups: `"all"` variants and `"on_target"` variants
- `TMB = n_variants / exome_target_size_Mb`
- ggplot2 bar plot with `geom_hline()` at 1, 10, 100 mutations/Mb (reference thresholds)

**Output:** `Tables/TMB.csv`, `Results/Variant_calling/TMB_per_territories.pdf`

---

### 7. MAF_ANALYSIS → `max_MAF_per_territories.pdf`

**Script:** `scripts/maf_plots.R`

Computes maximum allele frequency (FREQ) for four variant subsets per sample:
1. **all**: all variants
2. **on_target**: variants with `territory == "on_target"`
3. **selected**: variants with `selected == TRUE`
4. **selected_on_target**: selected AND on-target

Plots a ggplot2 point/line chart faceted by subset type.

**Output:** `Results/Variant_calling/max_MAF_per_territories.pdf`, `max_MAF_summary.csv`

---

### 8. CLINICAL_REPORTS → `<SAMPLE>_freq*_DP*.drug.xlsx`

**Script:** `scripts/clinical_reports_xlsx.R`

Produces one styled Excel workbook per sample:

1. Subsets to `selected == TRUE` variants (when `selected` column present)
2. Applies `FREQ >= min_freq` and `DP >= min_dp` thresholds
3. Combines `snp138` + `cosmic70` into a single `snp138_cosmic70` display column
4. Formats `FREQ` as `FREQ_pct` (percentage string)
5. Rounds `CADD13_PHRED` to 2 decimal places
6. Restricts columns to OUTCOL (see `clinical_reports_xlsx.R`)
7. Writes `.xlsx` with **bold header**, **blue fill**, **cell borders**, **auto-width columns** (using `openxlsx`)

**OUTCOL definition:**
`Chr`, `Start`, `End`, `Ref`, `Alt`, `Gene.refGene`, `Func.refGene`, `ExonicFunc.refGene`, `AAChange.refGene`, `snp138_cosmic70`, `FREQ_pct`, `CADD13_PHRED`, `DP`, `variant_type`, `mc`, `druggable`, `drug_interactions`, `territory`, `sample`

**Customisation:**
- Adjust thresholds: `--min_freq` (default 0.05) and `--min_dp` (default 10)
- Modify `OUTCOL` in `clinical_reports_xlsx.R` to add or remove columns

---

## Column schemas at each stage

| Stage | Key new columns added |
|-------|-----------------------|
| PARSE_PRETTY_CSV | `n_methods`, `DP.A`, `sample` |
| ANNOTATE_TARGET_BAIT | `on_target`, `on_bait`, `territory` |
| SELECT_CANDIDATES | `selected` |
| ANNOTATE_FUNCTIONAL_DRUG | `variant_type`, `mc`, `druggable`, `drug_interactions` |
| CLINICAL_REPORTS | `snp138_cosmic70`, `FREQ_pct` (display only) |
