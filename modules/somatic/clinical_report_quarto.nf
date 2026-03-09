// Render a Quarto (.qmd) clinical report from pipeline outputs.
// Input: somatic_mutations.rds + downstream QC/report outputs (optional MultiQC HTML)
// Output: clinical_report.html (self-contained HTML)

process CLINICAL_REPORT_QUARTO {
  tag "quarto_report"
  publishDir "${params.outdir}/Clinical_Reports/", mode: params.publish_mode
  input:
    path(somatic_rds)
    path(qc_pdf)
    path(tmb_csv)
    path(maf_pdf)
    path(clinical_csv)
    path(multiqc_html)
  output:
    path("clinical_report.html"), emit: report
  script:
    def multiqc_arg = (multiqc_html.name != "null" && multiqc_html.name != "empty.txt") ? "multiqc_html=${multiqc_html}" : ""
    """
    # Copy Quarto template and assets into the work directory
    cp ${projectDir}/docs/clinical_report_quarto.qmd  ./clinical_report.qmd
    cp ${projectDir}/docs/clinical-report.css         ./clinical-report.css   || true
    cp ${projectDir}/docs/clinical-report-extra.html  ./clinical-report-extra.html || true

    ${params.quarto} render clinical_report.qmd \\
      --output clinical_report.html \\
      --execute-params "somatic_rds=${somatic_rds},qc_pdf=${qc_pdf},tmb_csv=${tmb_csv},maf_pdf=${maf_pdf},clinical_csv=${clinical_csv},${multiqc_arg}"
    """
}
