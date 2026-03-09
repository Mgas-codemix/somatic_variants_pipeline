// BWA mem alignment.
// Input 1: (sample_id, r1, r2, patient_id, kit, project) from TRIM_FASTQ.
// Input 2 (optional): (ref_fa, *.amb, *.ann, *.bwt, *.pac, *.sa) from BUILD_BWA_INDEX; when provided, ref is taken from here.
// Output: SAM then converted to BAM (sort + index).

process ALIGN_BWA {
  tag "${sample_id}"
  input:
    tuple val(sample_id), path(r1), path(r2), val(patient_id), val(kit), val(project)
    tuple path(ref_fa), path(ref_amb), path(ref_ann), path(ref_bwt), path(ref_pac), path(ref_sa)
  output:
    tuple val(sample_id), path("*.bam"), path("*.bam.bai"), val(patient_id), val(kit), val(project), emit: bam
  script:
    def ref = ref_fa.name
    def pe = r2.name != "null" && r2.name != "empty.txt" && r2.size() > 0
    // Use \\t so script gets \t and bash expands to tab (must use double quotes around -R value)
    def readgroup = "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tLB:lib1\\tPL:ILLUMINA\\tPU:${kit}"
    if (pe) {
      """
      ${params.bwa} mem -t ${task.cpus} -R "${readgroup}" ${ref} ${r1} ${r2} | ${params.samtools} sort -@ ${task.cpus} -O BAM -o ${sample_id}.sorted.bam -
      ${params.samtools} index ${sample_id}.sorted.bam
      """
    } else {
      """
      ${params.bwa} mem -t ${task.cpus} -R "${readgroup}" ${ref} ${r1} | ${params.samtools} sort -@ ${task.cpus} -O BAM -o ${sample_id}.sorted.bam -
      ${params.samtools} index ${sample_id}.sorted.bam
      """
    }
}
