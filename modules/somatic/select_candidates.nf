// selectionOfCandidates + checkColnames; filter and standardize columns
// Input: annotated table; ref Rdata paths (falseCancerGenes, gene_categories_repo)
// Output: candidate table (RDS)

process SELECT_CANDIDATES {
  tag "candidates"
  publishDir "${params.outdir}/Rdata", mode: params.publish_mode, pattern: "*.rds"
  input:
    path(annotated_rds)
    val(ref_rdata_dir)
    val(false_cancer_genes)
    val(gene_categories_repo)
  output:
    path("candidates.rds"), emit: out
  script:
    def ref_dir = ref_rdata_dir ?: "."
    def fcg = false_cancer_genes ?: "${ref_dir}/falseCancerGenes.Rdata"
    def gcr = gene_categories_repo ?: "${ref_dir}/gene_categories_repo.Rdata"
    """
    Rscript ${projectDir}/scripts/selection_candidates.R \\
      --input ${annotated_rds} \\
      --false_cancer_genes ${fcg} \\
      --gene_categories_repo ${gcr} \\
      --output candidates.rds
    """
}
