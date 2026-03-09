# Test data for somatic variants pipeline

## Contents

- **mock_tumor_R1.fastq.gz**, **mock_tumor_R2.fastq.gz** – Small mock paired-end FASTQ (5 reads, 50 bp).
- **samplesheet_fastq.csv** – Samplesheet with one sample `mock_tumor` and `run_mode=from_fastq`.
- **target.bed**, **bait.bed** – Minimal 3-column BEDs for downstream (parse → annotate → QC/reports).

## Run the pipeline locally

From the **pipeline root** (parent of `test_data/`):

```bash
nextflow run main.nf -entry from_fastq \
  --samplesheet test_data/samplesheet_fastq.csv \
  --target_bed test_data/target.bed \
  --bait_bed test_data/bait.bed \
  --outdir results_test \
  -profile local
```

The **from_fastq** path needs reference data and tools. If they are not set, the run will stop at the first process that needs them (e.g. **TRIM_FASTQ** if cutadapt/seqtk are missing, or **ALIGN_BWA** if `ref_fasta` is not set).

For a full run you must also pass (e.g.):

- `--ref_fasta /path/to/genome.fa`
- `--bwa_index /path/to/genome.fa` (or same as ref_fasta)
- `--dbsnp_vcf`, `--cosmic_vcf`, `--annovar_basedir`
- `--gatk_jar`, `--picard_jar`, `--mutect_jar`, `--varscan_jar`

To only test that the workflow is wired and the samplesheet is read, the command above is enough; it will fail with a clear error when a required tool or parameter is missing.
