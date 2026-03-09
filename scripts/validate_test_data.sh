#!/usr/bin/env bash
# Validate mock FASTQ test data and samplesheet. Run from pipeline root.
set -e
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

echo "Validating test data..."
FAIL=0

for f in test_data/mock_tumor_R1.fastq.gz test_data/mock_tumor_R2.fastq.gz test_data/samplesheet_fastq.csv test_data/target.bed test_data/bait.bed; do
  if [ -f "$f" ]; then
    echo "  OK $f"
  else
    echo "  MISSING $f"
    FAIL=1
  fi
done

if [ ! -f assets/empty.txt ]; then
  echo "  MISSING assets/empty.txt"
  FAIL=1
else
  echo "  OK assets/empty.txt"
fi

if [ $FAIL -eq 1 ]; then
  echo "Validation failed. Fix missing files and re-run."
  exit 1
fi

echo "All test data present. Run the pipeline with:"
echo "  ./run_test.sh"
echo "Or:"
echo "  nextflow run main.nf -entry from_fastq \\"
echo "    --samplesheet test_data/samplesheet_fastq.csv \\"
echo "    --target_bed test_data/target.bed \\"
echo "    --bait_bed test_data/bait.bed \\"
echo "    --outdir results_test \\"
echo "    -profile local"
exit 0
