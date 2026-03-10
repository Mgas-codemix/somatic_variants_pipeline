// ExonicFunc refGene mapping; variant_type (SNV/InDel); merge DGIdb + actionable; druggable flag
// Output: somatic_mutations.rds (and optional TSV)
// Emits out_rds for downstream QC/reports

process ANNOTATE_FUNCTIONAL_DRUG {
  tag "functional_drug"
  publishDir "${params.outdir}/Rdata", mode: params.publish_mode, pattern: "*.rds"
  input:
    path(candidates_rds)
    val(acc_actionable)
    val(dgidb_tsv)
  output:
    path("somatic_mutations.rds"), emit: out
  script:
    def acc_path = acc_actionable ?: ""
    def dgidb_path = dgidb_tsv ?: ""
    """
    Rscript ${projectDir}/scripts/annotate_functional_drug.R \\
      --input=${candidates_rds} \\
      --acc_actionable="${acc_path}" \\
      --dgidb="${dgidb_path}" \\
      --output=somatic_mutations.rds
    """
}
