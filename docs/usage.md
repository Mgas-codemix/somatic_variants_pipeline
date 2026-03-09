# Usage Guide — Somatic Variants Pipeline

## Prerequisites

| Requirement | Version | Notes |
|-------------|---------|-------|
| Nextflow | ≥ 22.0 (DSL2) | `curl -s https://get.nextflow.io | bash` |
| Java | 11+ | Required by Nextflow |
| R | ≥ 4.0 | For downstream scripts |
| Docker (optional) | 20+ | Use `-profile docker` for fully reproducible runs |
| Conda (optional) | | Use `environment.yml` for a managed environment |

---

## Installation

### Option 1: Conda (recommended)

```bash
conda env create -f environment.yml
conda activate somatic-variants
```

This installs: Nextflow, Java, R (with ggplot2, ggsci, openxlsx, pROC, reshape2, data.table, GenomicRanges), Python, BWA, Samtools, Cutadapt, seqtk, FastQC.

For the full `from_fastq` pipeline you additionally need GATK 3.x, Picard, MuTect 1.x, VarScan, and ANNOVAR — install them separately or use Docker.

### Option 2: Docker (fully reproducible)

```bash
# Build the image (from repo root)
docker build -t somatic-variants:latest -f docker/Dockerfile .

# Run any pipeline command with -profile docker
nextflow run main.nf -profile docker -entry from_pretty_csv \
    --samplesheet samplesheet.example.csv \
    --target_bed /path/to/target.bed \
    --bait_bed   /path/to/bait.bed \
    --outdir     results
```

### Option 3: R packages only (downstream / from_pretty_csv)

```bash
Rscript scripts/install_r_packages.R
```

---

## Samplesheet format

All runs require `--samplesheet <path/to/samplesheet.csv>`.

The samplesheet is a CSV with a header row. Rows with `run_mode=from_pretty_csv` are used for downstream-only runs; rows with `run_mode=from_fastq` for full upstream runs.

### Column reference

| Column | Required | Description |
|--------|----------|-------------|
| `sample_id` | **Yes** | Unique sample identifier |
| `run_mode` | **Yes** | `from_pretty_csv` or `from_fastq` |
| `patient_id` | No | Patient/subject ID (defaults to `sample_id`) |
| `kit` | Yes | `WHOLE_EXOME` or `GENE_PANEL` |
| `project` | No | Project / cohort label |
| `variants_pretty_csv` | When `run_mode=from_pretty_csv` | Path to `*_variants_pretty.csv` |
| `fastq_r1` | When `run_mode=from_fastq` | Path to R1 FASTQ |
| `fastq_r2` | When `run_mode=from_fastq` | Path to R2 FASTQ (`.` for single-end) |
| `adapter_r1` / `adapter_r2` | No | Override Cutadapt adapters |
| `cancer_site` | No | e.g. `soft_tissue_sarcoma`, `lung` |

### Example — from_pretty_csv

```csv
sample_id,run_mode,kit,project,variants_pretty_csv
SAMPLE1,from_pretty_csv,WHOLE_EXOME,ProjectA,/data/SAMPLE1_variants_pretty.csv
SAMPLE2,from_pretty_csv,WHOLE_EXOME,ProjectA,/data/SAMPLE2_variants_pretty.csv
```

### Example — from_fastq

```csv
sample_id,run_mode,kit,patient_id,fastq_r1,fastq_r2
TUMOR01,from_fastq,WHOLE_EXOME,PATIENT01,/data/TUMOR01_R1.fastq.gz,/data/TUMOR01_R2.fastq.gz
```

---

## Running the pipeline

### Downstream only (from pretty CSV)

```bash
nextflow run main.nf \
  -entry from_pretty_csv \
  --samplesheet samplesheet.example.csv \
  --target_bed /path/to/target.bed \
  --bait_bed   /path/to/bait.bed \
  --outdir     results
```

### Full pipeline (from FASTQ)

```bash
nextflow run main.nf \
  -entry from_fastq \
  --samplesheet samplesheet.example.csv \
  --ref_fasta   /path/to/genome.fa \
  --dbsnp_vcf   /path/to/dbsnp.vcf.gz \
  --cosmic_vcf  /path/to/cosmic.vcf.gz \
  --annovar_basedir /path/to/annovar \
  --gatk_jar    /path/to/GenomeAnalysisTK.jar \
  --picard_jar  /path/to/picard.jar \
  --mutect_jar  /path/to/mutect-1.1.7.jar \
  --varscan_jar /path/to/VarScan.v2.4.4.jar \
  --target_bed  /path/to/target.bed \
  --bait_bed    /path/to/bait.bed \
  --outdir      results
```

### With optional benchmarking

Provide a truth VCF (e.g. GIAB) to run AUROC / PR curve evaluation:

```bash
nextflow run main.nf \
  -entry from_pretty_csv \
  --samplesheet samplesheet.example.csv \
  --target_bed  /path/to/target.bed \
  --bait_bed    /path/to/bait.bed \
  --truth_vcf   /path/to/HG001_GRCh38_benchmark.vcf.gz \
  --confident_bed /path/to/HG001_GRCh38_confident.bed \
  --run_benchmark \
  --outdir      results
```

---

## Parameter reference

### Input

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--samplesheet` | — | **Required.** Path to sample sheet CSV |
| `--variants_dir` | — | Alternative: directory of pretty CSVs |

### Reference data (downstream)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--target_bed` | — | Capture target BED (required) |
| `--bait_bed` | — | Bait / padded BED (required) |
| `--ref_rdata_dir` | — | Directory containing reference Rdata files |
| `--false_cancer_genes` | — | `falseCancerGenes.Rdata` path |
| `--gene_categories_repo` | — | `gene_categories_repo.Rdata` path |
| `--acc_actionable_164` | — | `ACC_actionable_164.Rdata` path |
| `--dgidb_tsv` | — | `DGIdb_interactions.tsv` path |

### Reference data (upstream / from_fastq)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--ref_fasta` | — | Reference genome FASTA |
| `--dbsnp_vcf` | — | dbSNP VCF (for BQSR and MuTect) |
| `--cosmic_vcf` | — | COSMIC VCF (for MuTect) |
| `--annovar_basedir` | — | ANNOVAR installation directory |
| `--annovar_buildver` | `hg38` | Genome build for ANNOVAR |

### Tool paths

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--bwa` | `bwa` | BWA executable (or full path) |
| `--samtools` | `samtools` | Samtools executable |
| `--java` | `java` | Java executable |
| `--gatk_jar` | — | GATK 3 JAR |
| `--picard_jar` | — | Picard JAR |
| `--mutect_jar` | — | MuTect JAR |
| `--varscan_jar` | — | VarScan JAR |
| `--cutadapt` | `cutadapt` | Cutadapt executable |
| `--fastqc` | `fastqc` | FastQC executable |
| `--seqtk` | `seqtk` | seqtk executable |

### Output and thresholds

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--outdir` | `./results` | Output directory |
| `--min_freq` | `0.05` | Min allele frequency for clinical reports |
| `--min_dp` | `10` | Min read depth for clinical reports |
| `--min_allele_frequency` | `0.01` | Min allele frequency in caller merge step |
| `--min_cov_position` | `10` | Min coverage in caller merge step |

### Benchmarking (opt-in)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--truth_vcf` | — | Truth set VCF (e.g. GIAB) |
| `--confident_bed` | — | Confident regions BED |
| `--run_benchmark` | `false` | Activate benchmarking |

---

## Troubleshooting FAQ

**Q: `Error: --samplesheet file not found`**  
A: Check that the path is correct and readable. Use absolute paths.

**Q: `Unable to access jarfile /opt/somatic/jars/picard.jar`**  
A: You are using `-profile docker` but the container is not built. Run:
```bash
docker build -t somatic-variants:latest -f docker/Dockerfile .
```

**Q: `R package 'GenomicRanges' not found`**  
A: Install via Conda (`bioconductor-genomicranges`) or run:
```r
BiocManager::install("GenomicRanges")
```
The pipeline falls back to data.table interval arithmetic if GenomicRanges is unavailable.

**Q: Outputs are in the wrong location**  
A: Check `--outdir` and the publishDir settings in `conf/modules.config`.

**Q: How do I resume a failed run?**  
A: Add `-resume` to the Nextflow command. Nextflow caches completed tasks.
