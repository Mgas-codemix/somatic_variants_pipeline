// MuTect 1.x tumor-only (or tumor-normal). Input: (sample_id, bam, bai, patient_id, kit, project)
// Output: MuTect VCF

process MUTECT {
  tag "${sample_id}"
  input:
    tuple val(sample_id), path(bam), path(bai), val(patient_id), val(kit), val(project)
    path ref_fasta
    path dbsnp_vcf
    path cosmic_vcf
  output:
    tuple val(sample_id), path("*.vcf"), val(patient_id), val(kit), val(project), emit: vcf
  script:
    """
    ${params.java} -Xmx4g -jar ${params.mutect_jar} --reference_sequence ${ref_fasta} --dbsnp ${dbsnp_vcf} --cosmic ${cosmic_vcf} --input_file:tumor ${bam} --out ${sample_id}.mutect.vcf
    """
}
