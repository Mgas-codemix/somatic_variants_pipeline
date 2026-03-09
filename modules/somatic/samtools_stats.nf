// Post-alignment QC: samtools stats/flagstat/idxstats on recalibrated BAM.
// Input: tuple (sample_id, bam, bai, patient_id, kit, project)
// Output: stats files for MultiQC

process SAMTOOLS_STATS {
  tag "${sample_id}"
  publishDir "${params.outdir}/SamtoolsStats/", mode: params.publish_mode
  input:
    tuple val(sample_id), path(bam), path(bai), val(patient_id), val(kit), val(project)
  output:
    path("*.stats"),    emit: stats
    path("*.flagstat"), emit: flagstat
    path("*.idxstats"), emit: idxstats
  script:
    """
    ${params.samtools} stats    ${bam} > ${sample_id}.stats
    ${params.samtools} flagstat ${bam} > ${sample_id}.flagstat
    ${params.samtools} idxstats ${bam} > ${sample_id}.idxstats
    """
}
