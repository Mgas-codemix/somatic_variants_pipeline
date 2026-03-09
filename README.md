# Somatic Variants Pipeline

Modular Nextflow pipeline for somatic variant processing and annotation: from FASTQ (or from pre-called pretty CSV) through alignment, calling, annotation, QC, and clinical reports.

## Run modes

| Entry point | Description |
|-------------|-------------|
| **from_fastq** | FASTQ → trim → align (BWA) → BAM prep → mark dup → BQSR → MuTect + VarScan → merge → ANNOVAR → pretty CSV → downstream (parse, annotate, QC, reports). |
| **from_pretty_csv** | Downstream only: parse existing `*_variants_pretty.csv` → annotate target/bait → candidates → functional/drug → QC (territories, TMB, MAF) → clinical Excel reports. |

## Samplesheet (CSV)

The pipeline is driven by a **samplesheet** CSV. All runs require `--samplesheet <path>`.

### Columns (see `schema/samplesheet_schema.csv`)

| Column | Required | Description |
|--------|----------|-------------|
| sample_id | Yes | Unique sample identifier |
| run_mode | Yes | `from_pretty_csv` or `from_fastq` |
| patient_id | No | Patient/subject ID for grouping |
| kit | Yes | `WHOLE_EXOME` or `GENE_PANEL` |
| project | No | Project/cohort label |
| variants_pretty_csv | When run_mode=from_pretty_csv | Path to `*_variants_pretty.csv` |
| fastq_r1 | When run_mode=from_fastq | Path to R1 FASTQ (or FASTQ.GZ) |
| fastq_r2 | When run_mode=from_fastq | Path to R2 FASTQ, or `.` / empty for single-end |
| bam | When run_mode=from_bam | Path to BAM (future use) |
| tumor_sample_id | No | Tumor sample_id for tumor–normal pairing |
| normal_sample_id | No | Normal sample_id for tumor–normal pairing |
| adapter_r1 / adapter_r2 | No | Override adapters for Cutadapt |
| cancer_site | No | Cancer site for filtering (e.g. soft_tissue_sarcoma, lung) |
| read_group_* | No | Override read group fields |

### Example samplesheets

- **from_pretty_csv**: `samplesheet.example.csv` (last row).
- **from_fastq**: `samplesheet.example.csv` (first rows) or use the same file with mixed run_mode rows.

You can mix rows: some rows `from_fastq`, some `from_pretty_csv`; the pipeline runs the workflow corresponding to the entry point you choose (see below).

## Installing required packages

### Option 1: Conda (recommended)

Install Nextflow, Java, R, and the alignment/QC tools in one environment:

```bash
conda env create -f environment.yml
conda activate somatic-variants
```

This provides: Nextflow, OpenJDK 8, R (with data.table, xlsx), Python 3, BWA, Samtools, Cutadapt, seqtk, FastQC. For GATK, Picard, MuTect, VarScan and ANNOVAR (needed for full **from_fastq**), install them separately or use the Docker image (`-profile docker`).

### Option 2: R packages only (downstream / from_pretty_csv)

If you already have R and only need the R dependencies:

```bash
Rscript scripts/install_r_packages.R
```

### Option 3: Docker

Use the pipeline image so all tools run inside the container (see [Docker (reproducible runs)](#docker-reproducible-runs) below).

## Quick start

### Downstream only (from pretty CSV)

```bash
cd somatic_variants_pipeline
nextflow run main.nf \
  -entry from_pretty_csv \
  --samplesheet samplesheet.example.csv \
  --target_bed /path/to/target.bed \
  --bait_bed /path/to/bait.bed \
  --outdir results
```

### Full pipeline (from FASTQ)

```bash
cd somatic_variants_pipeline
nextflow run main.nf \
  -entry from_fastq \
  --samplesheet samplesheet.example.csv \
  --ref_fasta /path/to/genome.fa \
  --bwa_index /path/to/genome.fa \
  --dbsnp_vcf /path/to/dbsnp.vcf.gz \
  --cosmic_vcf /path/to/cosmic.vcf.gz \
  --annovar_basedir /path/to/annovar \
  --gatk_jar /path/to/GenomeAnalysisTK.jar \
  --picard_jar /path/to/picard.jar \
  --mutect_jar /path/to/mutect-1.1.7.jar \
  --varscan_jar /path/to/VarScan.v2.4.4.jar \
  --target_bed /path/to/target.bed \
  --bait_bed /path/to/bait.bed \
  --outdir results
```

Set `--skip_cutadapt true` to skip adapter trimming and only run SeqTK trim.

## Parameters (main)

See `nextflow.config` and override as needed.

- **Input**: `samplesheet`, `variants_dir`
- **References (downstream)**: `target_bed`, `bait_bed`, `ref_rdata_dir`, `false_cancer_genes`, `gene_categories_repo`, `acc_actionable_164`, `dgidb_tsv`
- **References (upstream)**: `ref_fasta`, `bwa_index`, `dbsnp_vcf`, `cosmic_vcf`, `known_indels`, `annovar_basedir`, `annovar_buildver`, `adapter_r1`, `adapter_r2`
- **Tool paths**: `bwa`, `samtools`, `java`, `gatk_jar`, `picard_jar`, `mutect_jar`, `varscan_jar`, `cutadapt`, `fastqc`, `seqtk`
- **Output**: `outdir`, `clinical_reports_dir`, `tmb_table_path`
- **Thresholds**: `min_freq`, `min_dp`, `tmb_ref_lines`, `min_allele_frequency`, `min_cov_position`, `varscan_min_coverage_tumor`, `varscan_min_var_freq`
- **Flags**: `run_mode`, `skip_cutadapt`, `skip_fastqc`, `publish_mode`

## Output layout

- **from_fastq**: `results/Variants/WXS/<sample>_variants_pretty.csv`, then same downstream outputs as below.
- **from_pretty_csv** (or after from_fastq):  
  `results/Rdata/somatic_mutations.rds`,  
  `results/Results/Variant_calling/QCs/`,  
  `results/Results/Variant_calling/TMB_per_territories.pdf`,  
  `results/Tables/TMB.csv`,  
  `results/Results/Variant_calling/CLINICAL_REPORTS/*.xlsx`.

## Docker (reproducible runs)

A single Docker image provides all tools. Build and run with the **docker** profile:

```bash
# Build (from repo root)
docker build -t somatic-variants:latest -f docker/Dockerfile .

# Run with Docker: use -profile docker (tool paths are set automatically)
nextflow run main.nf -profile docker \
  -entry from_pretty_csv \
  --samplesheet samplesheet.example.csv \
  --target_bed /path/to/target.bed \
  --bait_bed /path/to/bait.bed \
  --outdir results
```

Reference files must be visible inside the container (e.g. mount a host directory and pass paths like `--ref_fasta /data/refs/genome.fa`). See **docker/README.md** for:

- What is included (Picard, GATK 3.8, VarScan, BWA, Samtools, Cutadapt, seqtk, R, etc.)
- Optional MuTect and ANNOVAR (license or registration)
- Mounting reference data and `docker.runOptions`

## Requirements

- Nextflow 22+ (DSL2)
- For **from_pretty_csv**: R (data.table, etc.), reference Rdata and BEDs as in the Variants_Analysis document.
- For **from_fastq**: BWA, Samtools, GATK 3.x, Picard, MuTect 1.x, VarScan, ANNOVAR, Python 3, R; reference genome, DBSNP, COSMIC, target/bait BEDs.
- Or use **-profile docker** and the `somatic-variants` image (see Docker section above).

## Plan and schema

- Pipeline design: see plan in `.cursor/plans/` (somatic variants Nextflow pipeline).
- Samplesheet columns: `schema/samplesheet_schema.csv`.
- Example: `samplesheet.example.csv`.
