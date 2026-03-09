// Picard MarkDuplicates. Input: (sample_id, bam, bai, patient_id, kit, project)
// Output: dedup BAM + index + duplicate metrics (for MultiQC)

process MARK_DUPLICATES {
  tag "${sample_id}"
  input:
    tuple val(sample_id), path(bam), path(bai), val(patient_id), val(kit), val(project)
  output:
    tuple val(sample_id), path("*.md.bam"), path("*.md.bam.bai"), val(patient_id), val(kit), val(project), emit: bam
    path("*.md.metrics.txt"), emit: metrics
  script:
    """
    ${params.java} -Xmx4g -jar ${params.picard_jar} MarkDuplicates I=${bam} O=${sample_id}.md.bam M=${sample_id}.md.metrics.txt CREATE_INDEX=true
    ${params.samtools} index ${sample_id}.md.bam
    """
}
