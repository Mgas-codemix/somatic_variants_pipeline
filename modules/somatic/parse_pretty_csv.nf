// Parse and merge *_variants_pretty.csv; compute n_methods, DP.A, add sample column
// Input: tuple ( sample_list, csv_files ); sample_list = [ [sample_id, file, patient_id, kit, project], ... ]
// Output: single merged table (RDS)

process PARSE_PRETTY_CSV {
  tag "${sample_list.size()} samples"
  publishDir "${params.outdir}/Rdata", mode: params.publish_mode, pattern: "*.rds"
  input:
    tuple val(sample_list), path(csv_files)
  output:
    path("merged_variants.rds"), emit: out
  script:
    def csv_paths = csv_files instanceof List ? csv_files.join(",") : csv_files.toString()
    def sample_ids = sample_list.collect { it[0] }.join(",")
    """
    Rscript ${projectDir}/scripts/parse_and_enrich.R \\
      --csv_paths "${csv_paths}" \\
      --sample_ids "${sample_ids}" \\
      --output merged_variants.rds
    """
}
