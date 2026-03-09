/*
 * subworkflows/local/upstream.nf
 * FASTQ → pretty CSV upstream subworkflow.
 */

nextflow.enable.dsl = 2

include { TRIM_FASTQ }            from "../../modules/somatic/trim_fastq.nf"
include { BUILD_BWA_INDEX }       from "../../modules/somatic/build_bwa_index.nf"
include { ALIGN_BWA }             from "../../modules/somatic/align_bwa.nf"
include { ADD_READ_GROUPS }        from "../../modules/somatic/add_read_groups.nf"
include { MARK_DUPLICATES }       from "../../modules/somatic/mark_duplicates.nf"
include { BQSR }                  from "../../modules/somatic/bqsr.nf"
include { MUTECT }                from "../../modules/somatic/mutect.nf"
include { VARSCAN }               from "../../modules/somatic/varscan.nf"
include { MERGE_CALLERS_ANNOVAR } from "../../modules/somatic/merge_callers_annovar.nf"
include { COLLECT_PRETTY }        from "../../modules/somatic/collect_pretty.nf"

workflow UPSTREAM {
    take:
        ch_fastq_rows   // [sample_id, r1, r2, patient_id, kit, project, adapter_r1, adapter_r2]
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
