// Merge MuTect + VarScan VCFs, prepare ANNOVAR input, run ANNOVAR. Input: (sample_id, mutect_vcf, varscan_snp, varscan_indel, patient_id, kit, project)
// Output: multianno file path for collect_pretty

process MERGE_CALLERS_ANNOVAR {
  tag "${sample_id}"
  publishDir "${params.outdir}/Variants/WXS", mode: params.publish_mode, pattern: "*_variants_pretty.csv"
  input:
    tuple val(sample_id), path(mutect_vcf), path(varscan_snp), path(varscan_indel), val(patient_id), val(kit), val(project)
  output:
    tuple val(sample_id), path("*.multianno.txt"), path("*.variants.csv"), val(patient_id), val(kit), val(project), emit: annovar
  script:
    """
    python3 ${projectDir}/scripts/merge_callers_prepare_annovar.py \\
      --sample_id ${sample_id} \\
      --mutect ${mutect_vcf} \\
      --varscan_snp ${varscan_snp} \\
      --varscan_indel ${varscan_indel} \\
      --min_af ${params.min_allele_frequency} \\
      --min_dp ${params.min_cov_position} \\
      --out_variants ${sample_id}.variants.csv \\
      --out_annovar_input ${sample_id}.annovar_input
    perl ${params.annovar_basedir}/table_annovar.pl ${sample_id}.annovar_input ${params.annovar_basedir}/humandb/ -buildver ${params.annovar_buildver} -protocol refGene,snp138,cosmic70,clinvar_20160302,dbnsfp30a,cadd13 -operation g,f,f,f,f,f -nastring NA -remove -v
    """
}
