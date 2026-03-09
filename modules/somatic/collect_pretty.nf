// Write _variants_pretty.csv from multianno + variants table (HaTSPiL collect_annotated_variants logic).
// Input: (sample_id, multianno, variants_csv, patient_id, kit, project)
// Output: (sample_id, pretty_csv path, patient_id, kit, project) for downstream PARSE

process COLLECT_PRETTY {
  tag "${sample_id}"
  publishDir "${params.outdir}/Variants/WXS", mode: params.publish_mode, pattern: "*_variants_pretty.csv"
  input:
    tuple val(sample_id), path(multianno), path(variants_csv), val(patient_id), val(kit), val(project)
  output:
    tuple val(sample_id), path("*_variants_pretty.csv"), val(patient_id), val(kit), val(project), emit: pretty
  script:
    """
    Rscript ${projectDir}/scripts/collect_pretty.R \\
      --sample_id ${sample_id} \\
      --multianno ${multianno} \\
      --variants_csv ${variants_csv} \\
      --output ${sample_id}_variants_pretty.csv
    """
}
