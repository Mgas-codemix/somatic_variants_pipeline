#!/usr/bin/env bash
# Run somatic variants pipeline with mock FASTQ test data (from_fastq).
# Execute from pipeline root. Requires Nextflow.
#
# Option A - With Docker (recommended): start Docker Desktop, then:
#   docker build -t somatic-variants:latest -f docker/Dockerfile .
#   ./run_test.sh -profile docker
#
# Option B - Without Docker: provide JAR paths (Picard, GATK, MuTect, VarScan):
#   ./run_test.sh --picard_jar /path/to/picard.jar --gatk_jar /path/to/GenomeAnalysisTK.jar ...
#
# Mock reference (ref.fa, dbsnp.vcf, cosmic.vcf) is in test_data/reference/; BWA index is built by the pipeline.

set -e
cd "$(dirname "$0")"

# Validate test data first
bash scripts/validate_test_data.sh

if ! command -v nextflow &>/dev/null; then
  echo "Nextflow not found. Install from https://www.nextflow.io/ (requires Java)."
  exit 1
fi

# Use mock reference params; override with -profile docker or pass JAR paths
nextflow run main.nf -entry from_fastq \
  -c test_data/params_from_fastq_mock.config \
  -profile local \
  "$@"
