#!/usr/bin/env bash
# Build mock reference and BWA index for testing from_fastq.
# Run from project root: bash scripts/build_mock_reference.sh

set -euo pipefail
REFDIR="test_data/reference"
mkdir -p "$REFDIR"
REFFA="${REFDIR}/ref.fa"

# Generate small human-like reference: chr1 50000 bp, chr2 25000 bp (matches test_data target/bait BED)
if [[ ! -f "$REFFA" ]]; then
  python3 - <<'PY'
def seq(n, base="ACGT"):
    return (base * (n // len(base) + 1))[:n]
with open("test_data/reference/ref.fa", "w") as f:
    f.write(">chr1\n")
    for i in range(0, 50000, 80):
        f.write(seq(min(80, 50000 - i)) + "\n")
    f.write(">chr2\n")
    for i in range(0, 25000, 80):
        f.write(seq(min(80, 25000 - i)) + "\n")
PY
  echo "Created ${REFFA}"
else
  echo "Using existing ${REFFA}"
fi

# BWA index (required for alignment)
if command -v bwa &>/dev/null; then
  bwa index -a bwtsw "$REFFA"
  echo "BWA index built for ${REFFA}"
else
  echo "WARNING: bwa not in PATH; run: bwa index -a bwtsw ${REFFA}"
fi

# Optional: compress and index VCFs (requires bgzip + tabix)
if command -v bgzip &>/dev/null && command -v tabix &>/dev/null; then
  (cd "$REFDIR" && bgzip -cf dbsnp.vcf > dbsnp.vcf.gz && tabix -p vcf dbsnp.vcf.gz)
  (cd "$REFDIR" && bgzip -cf cosmic.vcf > cosmic.vcf.gz && tabix -p vcf cosmic.vcf.gz)
  echo "Compressed and indexed dbsnp.vcf.gz and cosmic.vcf.gz"
fi
