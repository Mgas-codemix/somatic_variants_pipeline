/*
 * lib/WorkflowSomatic.groovy
 * Somatic-variants-specific validation and helper utilities (nf-core style).
 */
class WorkflowSomatic {

    /** Validate the parsed samplesheet rows for the chosen entry point. */
    static void validateSamplesheet(rows, entry, log) {
        if (!rows) {
            log.error "No rows found in samplesheet for entry point: ${entry}"
            System.exit(1)
        }
        rows.each { row ->
            def sid = row.sample_id ?: row.get("sample_id")
            if (!sid) log.warn "Row missing sample_id: ${row}"
            if (entry == "from_pretty_csv") {
                def csv = row.variants_pretty_csv ?: row.get("variants_pretty_csv")
                if (!csv) log.warn "Row for sample ${sid} missing variants_pretty_csv"
                else if (!file(csv).exists()) log.warn "variants_pretty_csv not found: ${csv}"
            }
            if (entry == "from_fastq") {
                def r1 = row.fastq_r1 ?: row.get("fastq_r1")
                if (!r1) log.warn "Row for sample ${sid} missing fastq_r1"
            }
        }
    }

    /** Return a summary string of input counts. */
    static String summariseSamplesheet(rows) {
        def n_fastq  = rows.count { it.run_mode == "from_fastq" }
        def n_pretty = rows.count { it.run_mode == "from_pretty_csv" }
        return "Samplesheet: ${rows.size()} rows  [from_fastq: ${n_fastq}  from_pretty_csv: ${n_pretty}]"
    }
}
