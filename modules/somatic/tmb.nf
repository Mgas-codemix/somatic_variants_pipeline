// TMB = counts / exome_target_Mb; plot with ref lines; export TMB table
// Input: somatic_mutations.rds; target BED; outdir; tmb_table_path
// Output: TMB plot PDF and Tables/TMB.csv

process TMB {
  tag "tmb"
  publishDir "${params.outdir}/Results/Variant_calling", mode: params.publish_mode, pattern: "TMB*.pdf"
  publishDir dir: file(params.tmb_table_path).parent, mode: params.publish_mode, pattern: "TMB.csv"
  input:
    path(somatic_rds)
    path(target_bed)
    val(outdir)
    val(tmb_table_path)
  output:
    path("TMB.csv"), emit: out
  script:
    def tmb_dir = file(tmb_table_path).parent
    """
    mkdir -p ${tmb_dir}
    Rscript ${projectDir}/scripts/tmb.R \\
      --input=${somatic_rds} \\
      --target_bed=${target_bed} \\
      --outdir=${outdir} \\
      --tmb_csv=TMB.csv \\
      --ref_lines="${params.tmb_ref_lines.join(",")}"
    """
}
