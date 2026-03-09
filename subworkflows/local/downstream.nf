/*
 * subworkflows/local/downstream.nf
 * pretty CSV → annotated RDS → QC / TMB / MAF / clinical reports downstream subworkflow.
 */

nextflow.enable.dsl = 2

include { PARSE_PRETTY_CSV }         from "../../modules/somatic/parse_pretty_csv.nf"
include { ANNOTATE_TARGET_BAIT }     from "../../modules/somatic/annotate_target_bait.nf"
include { SELECT_CANDIDATES }       from "../../modules/somatic/select_candidates.nf"
include { ANNOTATE_FUNCTIONAL_DRUG } from "../../modules/somatic/annotate_functional_drug.nf"
include { QC_TERRITORIES }          from "../../modules/somatic/qc_territories.nf"
include { TMB }                     from "../../modules/somatic/tmb.nf"
include { MAF_ANALYSIS }            from "../../modules/somatic/maf_analysis.nf"
include { CLINICAL_REPORTS }        from "../../modules/somatic/clinical_reports.nf"

workflow DOWNSTREAM {
    take:
        ch_samplesheet_pretty   // collected rows (as produced by PARSE_PRETTY_CSV input)
        target_bed
        bait_bed

    main:
        PARSE_PRETTY_CSV(ch_samplesheet_pretty)
        ANNOTATE_TARGET_BAIT(PARSE_PRETTY_CSV.out, target_bed, bait_bed)
        SELECT_CANDIDATES(
            ANNOTATE_TARGET_BAIT.out,
            params.ref_rdata_dir     ?: ".",
            params.false_cancer_genes ?: "",
            params.gene_categories_repo ?: ""
        )
        ANNOTATE_FUNCTIONAL_DRUG(
            SELECT_CANDIDATES.out,
            params.acc_actionable_164 ?: "",
            params.dgidb_tsv          ?: ""
        )
        ch_somatic = ANNOTATE_FUNCTIONAL_DRUG.out
        QC_TERRITORIES(ch_somatic, target_bed, params.outdir)
        TMB(ch_somatic, target_bed, params.outdir, params.tmb_table_path)
        MAF_ANALYSIS(ch_somatic, params.outdir)
        CLINICAL_REPORTS(ch_somatic, params.outdir, params.clinical_reports_dir,
                         params.min_freq, params.min_dp)

    emit:
        somatic_rds = ch_somatic
        qc_pdf      = QC_TERRITORIES.out
        tmb_csv     = TMB.out
        reports     = CLINICAL_REPORTS.out
}
