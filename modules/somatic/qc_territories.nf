// Counts by territory (on_target, on_bait, off); bar plot variants_on_territories.pdf
// Input: somatic_mutations.rds; target BED; outdir
// Output: PDF and optional summary CSV

process QC_TERRITORIES {
  tag "qc_territories"
  publishDir "${params.outdir}/Results/Variant_calling/QCs", mode: params.publish_mode
  input:
    path(somatic_rds)
    path(target_bed)
    val(outdir)
  output:
    path("variants_on_territories.pdf"), emit: out
  script:
    """
    Rscript ${projectDir}/scripts/qc_territories.R \\
      --input=${somatic_rds} \\
      --target_bed=${target_bed} \\
      --outdir=${outdir} \\
      --pdf=variants_on_territories.pdf
    """
}
