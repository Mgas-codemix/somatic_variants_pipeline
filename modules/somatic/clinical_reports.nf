// Subset selected; OUTCOL; format FREQ, CADD, snp138-cosmic70; write styled Excel per sample
// Input: somatic_mutations.rds; outdir; clinical_reports_dir; min_freq; min_dp
// Output: CLINICAL_REPORTS/*.xlsx

process CLINICAL_REPORTS {
  tag "reports"
  publishDir "${params.clinical_reports_dir}", mode: params.publish_mode, pattern: "*.{csv,xlsx}"
  input:
    path(somatic_rds)
    val(outdir)
    val(clinical_reports_dir)
    val(min_freq)
    val(min_dp)
  output:
    path("*.csv"), emit: out
  script:
    """
    mkdir -p ${clinical_reports_dir}
    Rscript ${projectDir}/scripts/clinical_reports_xlsx.R \\
      --input ${somatic_rds} \\
      --outdir ${clinical_reports_dir} \\
      --min_freq ${min_freq} \\
      --min_dp ${min_dp}
    """
}
