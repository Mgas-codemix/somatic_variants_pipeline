// Aggregate all QC outputs from the pipeline into a single MultiQC report.
// Input: collected QC files from all stages
// Output: multiqc_report.html and multiqc_data/ directory

process MULTIQC {
  tag "multiqc"
  publishDir "${params.outdir}/MultiQC/", mode: params.publish_mode
  input:
    path(qc_files)
  output:
    path("multiqc_report.html"), emit: report
    path("multiqc_data/"),       emit: data
  script:
    def config_arg = params.multiqc_config ? "--config ${params.multiqc_config}" : ""
    """
    ${params.multiqc} . ${config_arg} --force --outdir .
    """
}
