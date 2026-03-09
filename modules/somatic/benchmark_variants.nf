/*
 * BENCHMARK_VARIANTS
 * Compare pipeline VCF calls against a truth set using rtg vcfeval or hap.py.
 * Produces TP/FP/FN counts, precision/recall, AUROC, and PR-curve PDFs.
 *
 * Inputs:
 *   tuple val(sample_id), path(calls_vcf)  — pipeline calls (VCF)
 *   path(truth_vcf)                         — truth set VCF
 *   path(confident_bed)                     — confident regions BED
 *   path(ref_sdf)                           — RTG SDF reference (optional; only for vcfeval)
 *   path(pretty_csv)                        — optional pretty CSV for score-based curves
 *
 * Outputs:
 *   bench_results  — directory with summary CSV, ROC PDF, PR PDF
 */

process BENCHMARK_VARIANTS {
    tag "${sample_id}"
    label "process_medium"

    publishDir "${params.outdir}/Benchmarking/${sample_id}", mode: params.publish_mode ?: "copy"

    input:
        tuple val(sample_id), path(calls_vcf)
        path truth_vcf
        path confident_bed
        path ref_sdf
        path pretty_csv

    output:
        tuple val(sample_id), path("benchmark_out/"), emit: bench_results

    script:
    def have_sdf = (ref_sdf && ref_sdf.name != "NO_FILE") ? "--ref_sdf ${ref_sdf}" : ""
    def have_bed = (confident_bed && confident_bed.name != "NO_FILE") ? confident_bed : ""
    def have_csv = (pretty_csv && pretty_csv.name != "NO_FILE") ? "--pretty_csv ${pretty_csv}" : ""
    """
    #!/usr/bin/env bash
    set -euo pipefail
    mkdir -p benchmark_out

    # ---- Try vcfeval (RTG Tools) ----
    if command -v rtg >/dev/null 2>&1 && [ -n "${have_sdf}" ]; then
        rtg vcfeval \\
            -b ${truth_vcf} \\
            -c ${calls_vcf} \\
            ${have_bed ? "--evaluation-regions ${have_bed}" : ""} \\
            ${have_sdf} \\
            -o benchmark_out/vcfeval \\
        && cp benchmark_out/vcfeval/summary.txt benchmark_out/summary.txt \\
        || echo "vcfeval failed or not available"

    # ---- Try hap.py (Illumina/GA4GH) ----
    elif command -v hap.py >/dev/null 2>&1; then
        hap.py ${truth_vcf} ${calls_vcf} \\
            ${have_bed ? "-f ${have_bed}" : ""} \\
            -o benchmark_out/happy \\
        && cp benchmark_out/happy.summary.csv benchmark_out/summary.csv \\
        || echo "hap.py failed or not available"

    else
        echo "Neither rtg vcfeval nor hap.py found; skipping variant comparison" \\
            > benchmark_out/summary.txt
    fi

    # ---- R: AUROC + PR curves ----
    Rscript ${projectDir}/scripts/benchmark_performance.R \\
        --bench_dir=benchmark_out \\
        ${have_csv} \\
        --outdir=benchmark_out \\
        --sample=${sample_id}
    """
}
