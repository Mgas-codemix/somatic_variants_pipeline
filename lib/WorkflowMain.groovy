/*
 * lib/WorkflowMain.groovy
 * Pipeline-level validation utilities (nf-core style).
 */
class WorkflowMain {

    /** Print a standard header banner to the log. */
    static void printHelp(workflow, params, log) {
        log.info """
        ============================================================
         Somatic Variants Pipeline
         Version : ${workflow.manifest.version ?: 'dev'}
         Nextflow : ${workflow.nextflow.version}
        ============================================================
        Usage:
          nextflow run main.nf -entry from_pretty_csv --samplesheet <CSV> ...
          nextflow run main.nf -entry from_fastq      --samplesheet <CSV> ...

        Mandatory parameters:
          --samplesheet   Path to sample sheet CSV
          --target_bed    Path to target capture BED
          --bait_bed      Path to bait BED

        Run nextflow run main.nf --help for full parameter list.
        ============================================================
        """.stripIndent()
    }

    /** Validate that mandatory file parameters point to existing files. */
    static void validateParams(params, log) {
        def errors = []

        if (!params.samplesheet) {
            errors << "--samplesheet is required"
        } else if (!file(params.samplesheet).exists()) {
            errors << "--samplesheet file not found: ${params.samplesheet}"
        }

        if (errors) {
            log.error "Validation errors:\n  " + errors.join("\n  ")
            System.exit(1)
        }
    }
}
