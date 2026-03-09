# Mock reference data for testing from_fastq

Use these files to run the pipeline from FASTQ without real human reference data.

## Contents

| File | Description |
|------|-------------|
| `ref.fa` | Small reference: chr1 (50 kb), chr2 (25 kb). Matches `test_data/target.bed` and `bait.bed`. |
| `dbsnp.vcf` | Minimal dbSNP-style VCF (chr1/chr2) for GATK BQSR and MuTect. |
| `cosmic.vcf` | Minimal COSMIC-style VCF for MuTect. |

## One-time setup

### 1. BWA index

The pipeline **builds the BWA index automatically** from `--ref_fasta` when running `from_fastq`. You do not need to run `bwa index` yourself.

(To (re)create only the mock `ref.fa` and optionally compress VCFs, run `bash scripts/build_mock_reference.sh` from the project root.)

### 2. Optional: compress and index VCFs

GATK and MuTect accept uncompressed `.vcf`; if you prefer `.vcf.gz` (e.g. for tabix):

```bash
cd test_data/reference
bgzip -c dbsnp.vcf > dbsnp.vcf.gz && tabix -p vcf dbsnp.vcf.gz
bgzip -c cosmic.vcf > cosmic.vcf.gz && tabix -p vcf cosmic.vcf.gz
```

Then use `--dbsnp_vcf test_data/reference/dbsnp.vcf.gz` and `--cosmic_vcf test_data/reference/cosmic.vcf.gz` in the run command below.

## Run from FASTQ with mock reference

You still need **Picard**, **GATK**, **MuTect**, and **VarScan** JARs (or use Docker).

**Option A – Direct run (real paths to JARs):**

From the project root:

```bash
nextflow run main.nf -entry from_fastq \
  --samplesheet test_data/samplesheet_fastq.csv \
  --ref_fasta "$(pwd)/test_data/reference/ref.fa" \
  --dbsnp_vcf "$(pwd)/test_data/reference/dbsnp.vcf" \
  --cosmic_vcf "$(pwd)/test_data/reference/cosmic.vcf" \
  --target_bed "$(pwd)/test_data/target.bed" \
  --bait_bed "$(pwd)/test_data/bait.bed" \
  --picard_jar /path/to/picard.jar \
  --gatk_jar /path/to/GenomeAnalysisTK.jar \
  --mutect_jar /path/to/muTect.jar \
  --varscan_jar /path/to/VarScan.jar
```

**Option B – Docker (JARs inside image):**

Build the image, then run with reference data mounted (see `docker/README.md`). Point `--ref_fasta`, `--dbsnp_vcf`, `--cosmic_vcf`, `--target_bed`, `--bait_bed` to the paths **inside the container** (e.g. `/data/refs/ref.fa` if you mount refs at `/data/refs`).

## ANNOVAR

The merge step (`MERGE_CALLERS_ANNOVAR`) requires ANNOVAR and a built `humandb` (e.g. refGene, snp138, cosmic70, etc.). This mock reference uses contigs `chr1`/`chr2` and build `hg38`; for a full test you must either install ANNOVAR and build humandb for hg38, or use a Docker image that includes ANNOVAR.
