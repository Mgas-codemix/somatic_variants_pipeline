# Docker image for Somatic Variants Pipeline

Single image with all tools required for both **from_fastq** and **from_pretty_csv** runs.

## Included tools

- **Java 8** (OpenJDK)
- **Alignment & QC**: BWA, Samtools, Cutadapt, seqtk, FastQC
- **Java tools**: Picard 2.27.5, GATK 3.8-1, VarScan 2.4.4
- **R** (with data.table, xlsx)
- **Python 3**, **Perl**
- **ANNOVAR**: directory prepared at `/opt/annovar` (you must add the scripts and optionally humandb; see below)

## Optional: MuTect 1.1.7

MuTect is not redistributable via public download. To include it in the image:

1. Obtain `muTect-1.1.7.jar` (Broad Institute or build from source).
2. Place it in `docker/jars/` (e.g. `docker/jars/muTect-1.1.7.jar`).
3. Rebuild the image.

If you do not add it, mount the JAR at run time and set `--mutect_jar` to the path inside the container.

## Optional: ANNOVAR

ANNOVAR requires registration: https://www.openbioinformatics.org/annovar/

1. Download `annovar.latest.tar.gz` and extract it.
2. Either:
   - **Option A**: Copy the extracted `annovar` directory into the image (custom Dockerfile step), or
   - **Option B**: Mount the ANNOVAR directory when running Nextflow so the container sees it at `/opt/annovar`.

The pipeline expects `annovar_basedir` to contain `table_annovar.pl` and (for full annotation) `humandb/` with built databases (e.g. hg38). Download databases with:

```bash
perl annotate_variation.pl -buildver hg38 -downdb refGene humandb/
# ... other -downdb as needed
```

## Build

From the **pipeline repository root**:

```bash
docker build -t somatic-variants:latest -f docker/Dockerfile .
```

To include MuTect, add the JAR to `docker/jars/` first, then build.

## Run the pipeline with Docker

Use the **docker** profile so every process runs inside the container:

```bash
nextflow run main.nf -profile docker \
  -entry from_fastq \
  --samplesheet test_data/samplesheet_fastq.csv \
  --target_bed test_data/target.bed \
  --bait_bed test_data/bait.bed \
  --ref_fasta /data/refs/genome.fa \
  --bwa_index /data/refs/genome.fa \
  --dbsnp_vcf /data/refs/dbsnp.vcf.gz \
  --cosmic_vcf /data/refs/cosmic.vcf.gz \
  --annovar_basedir /opt/annovar \
  --outdir results
```

**Reference data and paths**: Paths passed to `--ref_fasta`, `--dbsnp_vcf`, etc. must be valid **inside** the container. Either:

1. **Mount reference data** when running Nextflow (e.g. bind host refs to `/data/refs` and use `--ref_fasta /data/refs/genome.fa`), or  
2. **Copy refs into the image** in a custom Dockerfile and use those paths.

Nextflow mounts the pipeline directory and the current work directory into the container; any other host path must be mounted explicitly (e.g. `-v /path/to/refs:/data/refs`). How you pass that depends on your executor: with `-profile docker` and local executor, Nextflow runs `docker run` for each task; you can set `docker.runOptions = '-v /path/to/refs:/data/refs'` in `nextflow.config` or a custom config so that volume is added to every run.

### Example: mounting refs in Nextflow config

Create e.g. `conf/docker-refs.config`:

```groovy
docker {
  runOptions = '-v /path/on/host/to/refs:/data/refs'
}
```

Then run with:

```bash
nextflow run main.nf -profile docker -c conf/docker-refs.config ...
```

and use `--ref_fasta /data/refs/genome.fa`, etc.

## Reproducibility

- **Image tag**: Pin the image tag in `conf/docker.config` (e.g. `somatic-variants:20240304`) and use the same tag for production.
- **Reference versions**: Document the reference genome, DBSNP, COSMIC, and ANNOVAR DB versions you use.
- **Params**: Use a `-params-file params.json` (or similar) and keep it under version control so runs are reproducible.
