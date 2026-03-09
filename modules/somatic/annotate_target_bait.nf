// Annotate variants with on_target / on_bait / off using BED overlap (get_on_target)
// Input: merged table from PARSE_PRETTY_CSV; target BED; bait BED
// Output: table with territory tags

process ANNOTATE_TARGET_BAIT {
  tag "target_bait"
  publishDir "${params.outdir}/Rdata", mode: params.publish_mode, pattern: "*.rds"
  input:
    path(merged_rds)
    path(target_bed)
    path(bait_bed)
  output:
    path("annotated_territory.rds"), emit: out
  script:
    """
    Rscript ${projectDir}/scripts/annotate_target.R \\
      --input ${merged_rds} \\
      --target_bed ${target_bed} \\
      --bait_bed ${bait_bed} \\
      --output annotated_territory.rds
    """
}
