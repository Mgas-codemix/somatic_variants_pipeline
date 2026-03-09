// FastQC on raw or trimmed FASTQs.
// Input: tuple (sample_id, fastq_r1, fastq_r2)
// Output: FastQC zip files (for MultiQC) and HTML reports

process FASTQC {
  tag "${sample_id}"
  publishDir "${params.outdir}/FastQC/", mode: params.publish_mode
  input:
    tuple val(sample_id), path(fastq_r1), path(fastq_r2)
  output:
    path("*.zip"),  emit: zip
    path("*.html"), emit: html
  script:
    def r2_arg = (fastq_r2.name != "null" && fastq_r2.name != "empty.txt") ? "${fastq_r2}" : ""
    """
    ${params.fastqc} ${fastq_r1} ${r2_arg} --threads ${task.cpus} -o .
    """
}
