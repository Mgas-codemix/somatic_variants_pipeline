// GATK BaseRecalibrator + PrintReads (BQSR). Input: (sample_id, bam, bai, patient_id, kit, project)
// Output: recalibrated BAM + index

process BQSR {
  tag "${sample_id}"
  input:
    tuple val(sample_id), path(bam), path(bai), val(patient_id), val(kit), val(project)
    path ref_fasta
    path dbsnp_vcf
  output:
    tuple val(sample_id), path("*.recal.bam"), path("*.recal.bam.bai"), val(patient_id), val(kit), val(project), emit: bam
  script:
    def known_abs = params.known_indels ? (params.known_indels.toString().startsWith('/') ? params.known_indels : file(projectDir).resolve(params.known_indels.toString()).toAbsolutePath().toString()) : null
    def known = known_abs ? "-known ${known_abs}" : ""
    """
    ${params.java} -Xmx4g -jar ${params.gatk_jar} -T BaseRecalibrator -R ${ref_fasta} -I ${bam} -knownSites ${dbsnp_vcf} ${known} -o ${sample_id}.recal.table -nct ${task.cpus}
    ${params.java} -Xmx4g -jar ${params.gatk_jar} -T PrintReads -R ${ref_fasta} -I ${bam} -BQSR ${sample_id}.recal.table -o ${sample_id}.recal.bam -nct ${task.cpus}
    ${params.samtools} index ${sample_id}.recal.bam
    """
}
