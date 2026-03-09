/*
 * workflows/somatic_variants.nf
 * Main workflow logic, extracted from main.nf for nf-core modularity.
 */

nextflow.enable.dsl = 2

include { PARSE_PRETTY_CSV }             from "../modules/somatic/parse_pretty_csv.nf"
include { ANNOTATE_TARGET_BAIT }         from "../modules/somatic/annotate_target_bait.nf"
include { SELECT_CANDIDATES }           from "../modules/somatic/select_candidates.nf"
include { ANNOTATE_FUNCTIONAL_DRUG }     from "../modules/somatic/annotate_functional_drug.nf"
include { QC_TERRITORIES }              from "../modules/somatic/qc_territories.nf"
include { TMB }                         from "../modules/somatic/tmb.nf"
include { MAF_ANALYSIS }                from "../modules/somatic/maf_analysis.nf"
include { CLINICAL_REPORTS }            from "../modules/somatic/clinical_reports.nf"
include { TRIM_FASTQ }                  from "../modules/somatic/trim_fastq.nf"
include { BUILD_BWA_INDEX }             from "../modules/somatic/build_bwa_index.nf"
include { ALIGN_BWA }                   from "../modules/somatic/align_bwa.nf"
include { ADD_READ_GROUPS }              from "../modules/somatic/add_read_groups.nf"
include { MARK_DUPLICATES }             from "../modules/somatic/mark_duplicates.nf"
include { BQSR }                        from "../modules/somatic/bqsr.nf"
include { MUTECT }                      from "../modules/somatic/mutect.nf"
include { VARSCAN }                     from "../modules/somatic/varscan.nf"
include { MERGE_CALLERS_ANNOVAR }       from "../modules/somatic/merge_callers_annovar.nf"
include { COLLECT_PRETTY }              from "../modules/somatic/collect_pretty.nf"

// Optional benchmarking
include { BENCHMARK_VARIANTS }          from "../modules/somatic/benchmark_variants.nf"

workflow SOMATIC_VARIANTS_FROM_PRETTY {
    take:
        ch_samplesheet_pretty   // collected samplesheet rows for from_pretty_csv
        target_bed
        bait_bed

    main:
        PARSE_PRETTY_CSV(ch_samplesheet_pretty)
        ANNOTATE_TARGET_BAIT(PARSE_PRETTY_CSV.out, target_bed, bait_bed)
        SELECT_CANDIDATES(
            ANNOTATE_TARGET_BAIT.out,
            params.ref_rdata_dir ?: ".",
            params.false_cancer_genes ?: "",
            params.gene_categories_repo ?: ""
        )
        ANNOTATE_FUNCTIONAL_DRUG(
            SELECT_CANDIDATES.out,
            params.acc_actionable_164 ?: "",
            params.dgidb_tsv ?: ""
        )
        ch_somatic = ANNOTATE_FUNCTIONAL_DRUG.out
        QC_TERRITORIES(ch_somatic, target_bed, params.outdir)
        TMB(ch_somatic, target_bed, params.outdir, params.tmb_table_path)
        MAF_ANALYSIS(ch_somatic, params.outdir)
        CLINICAL_REPORTS(ch_somatic, params.outdir, params.clinical_reports_dir, params.min_freq, params.min_dp)

    emit:
        somatic_rds = ch_somatic
        qc_pdf      = QC_TERRITORIES.out
        tmb_csv     = TMB.out
        reports     = CLINICAL_REPORTS.out
}

workflow SOMATIC_VARIANTS_FROM_FASTQ {
    take:
        ch_fastq_rows
        ch_ref_fasta
        ch_dbsnp
        ch_cosmic

    main:
        BUILD_BWA_INDEX(ch_ref_fasta)
        TRIM_FASTQ(ch_fastq_rows)
        ALIGN_BWA(TRIM_FASTQ.out, BUILD_BWA_INDEX.out.ref_index)
        ADD_READ_GROUPS(ALIGN_BWA.out)
        MARK_DUPLICATES(ADD_READ_GROUPS.out)
        BQSR(MARK_DUPLICATES.out, ch_ref_fasta, ch_dbsnp)
        MUTECT(BQSR.out, ch_ref_fasta, ch_dbsnp, ch_cosmic)
        VARSCAN(BQSR.out, ch_ref_fasta)
        ch_merge = MUTECT.out.join(VARSCAN.out).map { m, v ->
            [ m[0], m[1], v[1], v[2], m[2], m[3], m[4] ]
        }
        MERGE_CALLERS_ANNOVAR(ch_merge)
        COLLECT_PRETTY(MERGE_CALLERS_ANNOVAR.out)

    emit:
        pretty_csvs = COLLECT_PRETTY.out
}
