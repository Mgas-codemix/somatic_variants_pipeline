// Max MAF overall and in target; plot max_MAF_per_territories.pdf
// Input: somatic_mutations.rds; outdir
// Output: PDF and optional MAF CSV

process MAF_ANALYSIS {
  tag "maf"
  publishDir "${outdir}/Results/Variant_calling", mode: params.publish_mode
  input:
    path(somatic_rds)
    val(outdir)
  output:
    path("max_MAF_per_territories.pdf"), emit: out_pdf
    path("max_MAF_summary.csv"), optional: true, emit: out_csv
  script:
    """
    Rscript ${projectDir}/scripts/maf_plots.R \\
      --input ${somatic_rds} \\
      --outdir ${outdir} \\
      --pdf max_MAF_per_territories.pdf
    """
}
