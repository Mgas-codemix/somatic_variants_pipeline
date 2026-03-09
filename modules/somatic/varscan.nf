// VarScan mpileup2snp + mpileup2indel (tumor-only). Input: (sample_id, bam, bai, patient_id, kit, project)
// Output: SNP VCF and Indel VCF

process VARSCAN {
  tag "${sample_id}"
  input:
    tuple val(sample_id), path(bam), path(bai), val(patient_id), val(kit), val(project)
    path ref_fasta
  output:
    tuple val(sample_id), path("*.snp.vcf"), path("*.indel.vcf"), val(patient_id), val(kit), val(project), emit: vcf
  script:
    """
    ${params.samtools} mpileup -d 8000 -f ${ref_fasta} ${bam} | ${params.java} -Xmx2g -jar ${params.varscan_jar} mpileup2snp --min-coverage ${params.varscan_min_coverage_tumor} --min-var-freq ${params.varscan_min_var_freq} --p-value 1 --output-vcf 1 > ${sample_id}.snp.vcf
    ${params.samtools} mpileup -d 8000 -f ${ref_fasta} ${bam} | ${params.java} -Xmx2g -jar ${params.varscan_jar} mpileup2indel --min-coverage ${params.varscan_min_coverage_tumor} --min-var-freq ${params.varscan_min_var_freq} --p-value 1 --output-vcf 1 > ${sample_id}.indel.vcf
    """
}
