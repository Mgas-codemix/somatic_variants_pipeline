# Variants_Analysis (1) – Pipeline mapping and spec

This document maps the **Clinical reports generation (WXS somatic)** section of `Variants_Analysis (1).docx` to the somatic variants Nextflow pipeline and lists outputs, parsing rules, and recommended enhancements.

---

## 1. Inputs (from doc)

| Source | Description |
|--------|-------------|
| **Variants/WXS** | CSV files matching pattern `*pretty*`, format `.csv` |
| **falseCancerGenes.Rdata** | From Pipelines/data (or ref_rdata_dir) |
| **gene_categories_repo.Rdata** | From Pipelines/data |
| **Target / bait BEDs** | e.g. Capture_kits/.../on_target.bed, Padded.bed |
| **ACC_actionable_164.Rdata** | Actionable genes list |
| **DGIdb_interactions.tsv** | Drug–gene interactions |

**Pipeline:** These map to `--target_bed`, `--bait_bed`, `--ref_rdata_dir` (or `--false_cancer_genes`, `--gene_categories_repo`), `--acc_actionable_164`, `--dgidb_tsv`, and `*_variants_pretty.csv` from COLLECT_PRETTY or provided for `from_pretty_csv`.

---

## 2. Parsing and enrichment (pretty CSV → merged table)

| Doc requirement | Implementation | Script / process |
|-----------------|----------------|------------------|
| **n_methods** | Count of methods = number of elements in `method` column (split by `:`) | `parse_and_enrich.R` (already: `lengths(strsplit(..., ":"))`) |
| **DP.A** | Allele-specific depth = FREQ × DP | `parse_and_enrich.R` (already) |
| **sample** | Add sample column per dataframe | `parse_and_enrich.R` (already) |
| Combine all sample data frames | `rbindlist(out_list, fill = TRUE)` | `parse_and_enrich.R` (already) |

**Output:** Single merged table → **merged_variants.rds** (PARSE_PRETTY_CSV).

---

## 3. Annotate target/bait and select candidates

| Doc requirement | Implementation | Script / process |
|-----------------|----------------|------------------|
| **get_on_target()** | Annotate variants as on_target / on_bait / off using BED overlap | `annotate_target.R` (MVP: sets territory; full: GenomicRanges/bedr overlap) |
| **selectionOfCandidates** | Filter variants using false cancer genes and gene categories | `selection_candidates.R` (MVP: loads Rdata, placeholder filter) |
| **checkColnames** | Ensure correct column names | `selection_candidates.R` (to be extended) |

**Output:** **annotated_territory.rds** → **candidates.rds**.

---

## 4. Functional and drug annotation

| Doc requirement | Implementation | Script / process |
|-----------------|----------------|------------------|
| **ExonicFunc.refGene** | Map to standard categories (mc) | `annotate_functional_drug.R` (MVP: variant_type only) |
| **variant_type** | SNV vs InDel from Ref/Alt length | `annotate_functional_drug.R` (already) |
| **DGIdb + actionable** | Merge drug interactions; aggregate by gene; match to variants | `annotate_functional_drug.R` (MVP: reads DGIdb, druggable = NA) |
| **druggable** | Flag from drug match / actionable criteria | `annotate_functional_drug.R` (to be filled) |
| Save **Rdata/somatic_mutations.rds** | Save final annotated table | `annotate_functional_drug.R` (already) |

**Output:** **somatic_mutations.rds** (and optional TSV).

---

## 5. QC – Variants by territory

| Doc requirement | Implementation | Script / process |
|-----------------|----------------|------------------|
| **ddply** | Counts by region (on_target, on_bait, off) per sample | `qc_territories.R` (already: `table(territory, sample)`) |
| **type** | Differentiate “all” vs “selected” variants | `qc_territories.R` (MVP: single type; doc: facet by type) |
| **ggplot** | Bar plot, facet_wrap(~type), scale_y_log10, ggsci/ggpubr | `qc_territories.R` (MVP: base barplot) |
| **Output path** | Results/Variant_calling/QCs/variants_on_territories.pdf | Process publishDir (already) |

**Enhancement:** Use `selected` column when present; two facets (all vs selected). Add optional ggplot2/ggsci/ggpubr for doc-style plots.

---

## 6. TMB (tumor mutational burden)

| Doc requirement | Implementation | Script / process |
|-----------------|----------------|------------------|
| **exome_target_size_Mb** | Total size of exome target regions from BED (Mb) | `tmb.R` (MVP: fixed 50 Mb; doc: sum BED intervals) |
| **stat** | Counts by region; TMB = counts / exome_target_size_Mb | `tmb.R` (MVP: total n / 50) |
| **geom_hline** | Reference lines at 1, 10, 100 | `tmb.R` (already) |
| **Tables/TMB.csv** | Export TMB table | `tmb.R` (already) |
| **Google Sheet** | Optional sheet_write (db_googlesheet.R) | Not in pipeline; optional config |

**Enhancement:** Compute exome_target_size_Mb from target BED (sum of interval lengths) and use it in TMB; optionally by territory (on_target, etc.).

---

## 7. MAF (max allele frequency)

| Doc requirement | Implementation | Script / process |
|-----------------|----------------|------------------|
| **maf** | Max allele frequency for “all” and “selected”, overall and in target | `maf_plots.R` (MVP: max FREQ per sample only) |
| **reshape2::melt** | Long format for plotting | `maf_plots.R` (to add) |
| **ggplot** | Line/point, facet_wrap(~type) | `maf_plots.R` (MVP: barplot) |
| **Output** | Results/Variant_calling/max_MAF_per_territories.pdf | Process publishDir (already) |

**Enhancement:** Compute max MAF (1) all variants, (2) selected only, (3) on-target only; facet by type; optional ggplot2.

---

## 8. Clinical reports (Excel)

| Doc requirement | Implementation | Script / process |
|-----------------|----------------|------------------|
| **readRDS("Rdata/somatic_mutations.rds")** | Load final annotated table | `clinical_reports_xlsx.R` (already) |
| **subset(df, selected)** | Only report selected variants | `clinical_reports_xlsx.R` (MVP: filter by FREQ/DP only; add selected when present) |
| **FREQ** | As percentage | `clinical_reports_xlsx.R` (format for display) |
| **snp138-cosmic70** | Single combined column | `clinical_reports_xlsx.R` (to add) |
| **OUTCOL** | Fixed set of columns for report | See below |
| **CADD13_PHRED** | Rounded for readability | `clinical_reports_xlsx.R` (to add) |
| **sample** | Clean ID (e.g. remove path/extension) | `clinical_reports_xlsx.R` (partial) |
| **createWorkbook / createSheet / CellStyle / addDataFrame / saveWorkbook** | Styled Excel | `clinical_reports_xlsx.R` (MVP: xlsx::write.xlsx; doc: bold headers, borders) |
| **Output** | Results/Variant_calling/CLINICAL_REPORTS/*_freq*_DP*.drug.xlsx | Process publishDir (already) |

**OUTCOL (from doc):** First 11 columns renamed to a standard set; typically include Chr, Start, End, Ref, Alt, Gene, Func, ExonicFunc, AAChange, snp138-cosmic70 (combined), FREQ, CADD, DP, and drug-related columns. Exact names can match existing pretty CSV + drug columns.

**Enhancement:** Define OUTCOL vector; subset to `selected == TRUE` when present; add snp138-cosmic70 column; round CADD13_PHRED; restrict written columns to OUTCOL; optional xlsx styling (bold header, borders).

---

## 9. Output layout (doc vs pipeline)

| Doc path | Pipeline (params.outdir) |
|----------|---------------------------|
| Variants/WXS/*_variants_pretty.csv | outdir/Variants/WXS/ |
| Rdata/somatic_mutations.rds | outdir/Rdata/ |
| Results/Variant_calling/QCs/variants_on_territories.pdf | outdir/Results/Variant_calling/QCs/ |
| Results/Variant_calling/TMB_per_territories.pdf | outdir/Results/Variant_calling/ |
| Results/Variant_calling/max_MAF_per_territories.pdf | outdir/Results/Variant_calling/ |
| Results/Variant_calling/CLINICAL_REPORTS/*.xlsx | params.clinical_reports_dir |
| Tables/TMB.csv | params.tmb_table_path (e.g. outdir/Tables/TMB.csv) |

Pipeline already matches this layout.

---

## 10. Summary of recommended script changes

1. **parse_and_enrich.R** – Already implements n_methods, DP.A, sample; no change required for doc.
2. **annotate_target.R** – Add real BED overlap (e.g. GenomicRanges or bedr) to set on_target / on_bait / off and territory.
3. **selection_candidates.R** – Implement selectionOfCandidates logic and checkColnames using loaded falseCancerGenes and gene_categories_repo; add/ensure `selected` column.
4. **annotate_functional_drug.R** – Add ExonicFunc mapping (mc); full DGIdb/actionable merge and druggable flag.
5. **qc_territories.R** – Optional: facet by type (all vs selected); optional ggplot2/ggsci/ggpubr and log scale.
6. **tmb.R** – Compute exome_target_size_Mb from target BED; use it in TMB; optionally TMB by territory.
7. **maf_plots.R** – Max MAF for all vs selected and overall vs on-target; facet by type; optional ggplot2.
8. **clinical_reports_xlsx.R** – Subset by `selected` when present; define OUTCOL; add snp138-cosmic70 column; round CADD13_PHRED; restrict to OUTCOL; optional xlsx styling (bold, borders).

This keeps the pipeline aligned with Variants_Analysis while leaving optional (e.g. Google Sheet, full ggplot2 themes) as config or future work.
