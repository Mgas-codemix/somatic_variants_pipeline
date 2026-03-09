// Picard AddOrReplaceReadGroups. Input: (sample_id, bam, bai, patient_id, kit, project)
// Output: RG BAM + index (used for mark_dup and downstream)

process ADD_READ_GROUPS {
  tag "${sample_id}"
  input:
    tuple val(sample_id), path(bam), path(bai), val(patient_id), val(kit), val(project)
  output:
    tuple val(sample_id), path("*.rg.bam"), path("*.rg.bam.bai"), val(patient_id), val(kit), val(project), emit: bam
  script:
    """
    ${params.java} -Xmx2g -jar ${params.picard_jar} AddOrReplaceReadGroups I=${bam} O=${sample_id}.rg.bam RGID=${sample_id} RGLB=lib1 RGPL=ILLUMINA RGSM=${sample_id} RGPU=${kit} CREATE_INDEX=true
    ${params.samtools} index ${sample_id}.rg.bam
    """
}
