#!/usr/bin/env bash
# Validate that all expected pipeline outputs exist and are non-empty after a
# from_pretty_csv pipeline run.
#
# Usage:
#   bash scripts/validate_test_outputs.sh [outdir]
#
# Default outdir: results_test_pretty

set -euo pipefail
OUTDIR="${1:-results_test_pretty}"

echo "=== Validating pipeline outputs in: ${OUTDIR} ==="
FAIL=0

check_exists() {
  local pattern="$1"
  local desc="$2"
  # Use mapfile to capture glob matches safely (handles spaces in filenames)
  local -a matches
  # shellcheck disable=SC2206
  mapfile -d '' matches < <(compgen -G "${OUTDIR}/${pattern}" 2>/dev/null | tr '\n' '\0' || true)
  if [ ${#matches[@]} -eq 0 ]; then
    echo "  MISSING  : ${desc} (${OUTDIR}/${pattern})"
    FAIL=1
  else
    for f in "${matches[@]}"; do
      if [ -s "$f" ]; then
        echo "  OK       : $f"
      else
        echo "  EMPTY    : $f"
        FAIL=1
      fi
    done
  fi
}

# --- Rdata / RDS outputs ---
check_exists "Rdata/merged_variants.rds"          "merged variants RDS"
check_exists "Rdata/somatic_mutations.rds"        "somatic mutations RDS"

# --- TMB outputs ---
# TMB.csv is published to <outdir>/Tables/ (params.tmb_table_path parent)
check_exists "Tables/TMB.csv"                      "TMB CSV"
# TMB PDF is published to <outdir>/Results/Variant_calling/
check_exists "Results/Variant_calling/TMB*.pdf"    "TMB PDF"

# --- QC territories ---
# variants_on_territories.pdf is published to <outdir>/Results/Variant_calling/QCs/
check_exists "Results/Variant_calling/QCs/variants_on_territories.pdf" "QC territories PDF"

# --- MAF analysis ---
# max_MAF_per_territories.pdf is published to <outdir>/Results/Variant_calling/
check_exists "Results/Variant_calling/max_MAF_per_territories.pdf" "MAF analysis PDF"

# --- Clinical reports ---
# Clinical reports are published to <outdir>/Variant_calling/CLINICAL_REPORTS/
check_exists "Variant_calling/CLINICAL_REPORTS/*.csv" "clinical report CSVs"

# --- Exit ---
echo ""
if [ $FAIL -eq 0 ]; then
  echo "All expected outputs are present and non-empty."
  exit 0
else
  echo "ERROR: One or more expected outputs are missing or empty."
  echo "       Check the Nextflow log and ${OUTDIR}/ for details."
  exit 1
fi
