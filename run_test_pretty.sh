#!/usr/bin/env bash
# Run the somatic variants pipeline with synthetic test data (from_pretty_csv downstream path).
# Usage:  ./run_test_pretty.sh [extra nextflow args]
set -euo pipefail
cd "$(dirname "$0")"

# Generate test Rdata reference files if they don't exist yet
if [ ! -f test_data/falseCancerGenes.Rdata ] || \
   [ ! -f test_data/gene_categories_repo.Rdata ] || \
   [ ! -f test_data/ACC_actionable_164.Rdata ]; then
    echo "Generating test reference Rdata files..."
    Rscript test_data/generate_test_rdata.R
fi

if ! command -v nextflow &>/dev/null; then
    echo "ERROR: nextflow not found. Install from https://www.nextflow.io/" >&2
    exit 1
fi

nextflow run main.nf -entry from_pretty_csv \
    -c test_data/params_from_pretty_csv_test.config \
    -profile local \
    "$@"
