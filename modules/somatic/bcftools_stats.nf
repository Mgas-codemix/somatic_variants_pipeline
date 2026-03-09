// Post-variant-calling QC: bcftools stats on caller VCF.
// Input: tuple (sample_id, vcf)
// Output: stats file for MultiQC

process BCFTOOLS_STATS {
  tag "${sample_id}"
  publishDir "${params.outdir}/BcftoolsStats/", mode: params.publish_mode
  input:
    tuple val(sample_id), path(vcf)
  output:
    path("*.bcftools_stats.txt"), emit: stats
  script:
    """
    ${params.bcftools} stats ${vcf} > ${sample_id}.bcftools_stats.txt
    """
}
