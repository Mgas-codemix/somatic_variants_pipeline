/*
 * Somatic Variants Pipeline - Main workflow
 * Entry points:
 *   - from_fastq: FASTQ -> BAM -> callers -> ANN OVAR -> pretty CSV (HaTSPiL-aligned upstream, stops at pretty.csv)
 *   - from_pretty_csv: pretty CSV -> downstream (Variants_Analysis-style QC, TMB, MAF, clinical reports)
 */

nextflow.enable.dsl = 2

def samplesheet = params.samplesheet
def target_bed  = params.target_bed
def bait_bed    = params.bait_bed
def outdir      = params.outdir

if (!samplesheet || !file(samplesheet).exists()) {
  exit 1, "ERROR: Provide a valid samplesheet CSV (--samplesheet <path>). See samplesheet.example.csv and schema/samplesheet_schema.csv"
}

// Validate from_fastq requirements so we fail fast with a clear message (skip when using container paths, e.g. -profile docker)
def hasFromFastq = file(samplesheet).readLines().drop(1).any { it.contains('from_fastq') }
def useContainerPaths = params.picard_jar?.toString()?.startsWith('/opt/')
if (hasFromFastq && useContainerPaths) {
  // Docker profile: ensure the container exists and has JARs so we don't fail later with "Unable to access jarfile"
  def dockerCheck = ["docker", "run", "--rm", "somatic-variants:latest", "test", "-f", "/opt/somatic/jars/picard.jar"].execute()
  dockerCheck.waitForOrKill(15000)
  if (dockerCheck.exitValue() != 0) {
    exit 1, """ERROR: from_fastq is using Docker profile but the somatic-variants container is missing or incomplete.
  (Processes run on the host with container paths, hence 'Unable to access jarfile'.)
  1. Start Docker Desktop and wait until it is running.
  2. Build the image: docker build -t somatic-variants:latest -f docker/Dockerfile .
  3. Re-run: nextflow run main.nf -entry from_fastq -profile docker -c test_data/params_from_fastq_mock.config -resume"""
  }
}
if (hasFromFastq && !useContainerPaths) {
  def missing = []
  def pj = params.picard_jar?.toString()?.replaceAll(/^["']|["']$/, '')
  if (!pj || pj == 'null' || pj.contains('path/to') || !file(pj).exists()) missing << "--picard_jar (path to picard.jar; use a real path or -profile docker)"
  def gj = params.gatk_jar?.toString()?.replaceAll(/^["']|["']$/, '')
  if (!gj || gj == 'null' || gj.contains('path/to') || !file(gj).exists()) missing << "--gatk_jar (path to GenomeAnalysisTK.jar)"
  def rf = params.ref_fasta?.toString()?.replaceAll(/^["']|["']$/, '')
  if (!rf || rf == 'null' || rf.contains('path/to') || !file(rf).exists()) missing << "--ref_fasta (reference FASTA)"
  if (missing) exit 1, "ERROR: from_fastq requires these existing files:\n  " + missing.join("\n  ")
}

// --- Load modules ---
include { PARSE_PRETTY_CSV }             from "./modules/somatic/parse_pretty_csv.nf"
include { ANNOTATE_TARGET_BAIT }         from "./modules/somatic/annotate_target_bait.nf"
include { SELECT_CANDIDATES }           from "./modules/somatic/select_candidates.nf"
include { ANNOTATE_FUNCTIONAL_DRUG }     from "./modules/somatic/annotate_functional_drug.nf"
include { QC_TERRITORIES }              from "./modules/somatic/qc_territories.nf"
include { TMB }                         from "./modules/somatic/tmb.nf"
include { MAF_ANALYSIS }                from "./modules/somatic/maf_analysis.nf"
include { CLINICAL_REPORTS }            from "./modules/somatic/clinical_reports.nf"
include { TRIM_FASTQ }                  from "./modules/somatic/trim_fastq.nf"
include { BUILD_BWA_INDEX }             from "./modules/somatic/build_bwa_index.nf"
include { ALIGN_BWA }                   from "./modules/somatic/align_bwa.nf"
include { ADD_READ_GROUPS }              from "./modules/somatic/add_read_groups.nf"
include { MARK_DUPLICATES }             from "./modules/somatic/mark_duplicates.nf"
include { BQSR }                        from "./modules/somatic/bqsr.nf"
include { MUTECT }                      from "./modules/somatic/mutect.nf"
include { VARSCAN }                     from "./modules/somatic/varscan.nf"
include { MERGE_CALLERS_ANNOVAR }       from "./modules/somatic/merge_callers_annovar.nf"
include { COLLECT_PRETTY }              from "./modules/somatic/collect_pretty.nf"
include { FASTQC as FASTQC_RAW }        from "./modules/somatic/fastqc.nf"
include { FASTQC as FASTQC_TRIMMED }    from "./modules/somatic/fastqc.nf"
include { SAMTOOLS_STATS }              from "./modules/somatic/samtools_stats.nf"
include { BCFTOOLS_STATS as BCFTOOLS_STATS_MUTECT }  from "./modules/somatic/bcftools_stats.nf"
include { BCFTOOLS_STATS as BCFTOOLS_STATS_VARSCAN } from "./modules/somatic/bcftools_stats.nf"
include { MULTIQC }                     from "./modules/somatic/multiqc.nf"
include { CLINICAL_REPORT_QUARTO }      from "./modules/somatic/clinical_report_quarto.nf"

// --- Channels from samplesheet ---
ch_samplesheet_pretty = Channel
  .fromPath(samplesheet)
  .splitCsv(header: true, strip: true)
  .filter { it.run_mode == "from_pretty_csv" }
  .map { row -> tuple(row.sample_id, file(row.variants_pretty_csv?.trim()), row.patient_id ?: row.sample_id, row.kit ?: "WHOLE_EXOME", row.project ?: "default") }
  .toList()
  .map { rows -> tuple(rows, rows.collect { it[1] }) }

ch_fastq_rows = Channel
  .fromPath(samplesheet)
  .splitCsv(header: true, strip: true)
  .filter { it.run_mode == "from_fastq" }
  .map { row ->
    def r1 = row.fastq_r1?.trim()
    def r2 = row.fastq_r2?.trim()
    def r1file = file(r1.startsWith('/') ? r1 : "${projectDir}/${r1}")
    def r2file = (r2 != null && r2 != "." && r2 != "") ? file(r2.startsWith('/') ? r2 : "${projectDir}/${r2}") : file("${projectDir}/assets/empty.txt")
    [ row.sample_id,
      r1file,
      r2file,
      row.patient_id ?: row.sample_id,
      row.kit ?: "WHOLE_EXOME",
      row.project ?: "default",
      row.adapter_r1 ?: params.adapter_r1,
      row.adapter_r2 ?: params.adapter_r2
    ]
  }

// --- Downstream (shared): parse -> annotate -> QC -> reports ---
workflow downstream_from_parse {
  take:
    parse_out
    multiqc_html  // optional: path to MultiQC report (empty file when not available)
  main:
    ANNOTATE_TARGET_BAIT(parse_out, file(target_bed), file(bait_bed))
    SELECT_CANDIDATES(ANNOTATE_TARGET_BAIT.out, params.ref_rdata_dir ?: ".", params.false_cancer_genes ?: "", params.gene_categories_repo ?: "")
    ANNOTATE_FUNCTIONAL_DRUG(SELECT_CANDIDATES.out, params.acc_actionable_164 ?: "", params.dgidb_tsv ?: "")
    ch_somatic_rds = ANNOTATE_FUNCTIONAL_DRUG.out
    QC_TERRITORIES(ch_somatic_rds, file(target_bed), params.outdir)
    TMB(ch_somatic_rds, file(target_bed), params.outdir, params.tmb_table_path)
    MAF_ANALYSIS(ch_somatic_rds, params.outdir)
    CLINICAL_REPORTS(ch_somatic_rds, params.outdir, params.clinical_reports_dir, params.min_freq, params.min_dp)
    ch_qc_pdf   = QC_TERRITORIES.out
    ch_tmb_csv  = TMB.out
    ch_reports  = CLINICAL_REPORTS.out
    if (!params.skip_quarto_report) {
      CLINICAL_REPORT_QUARTO(
        ch_somatic_rds,
        ch_qc_pdf,
        ch_tmb_csv,
        MAF_ANALYSIS.out.out_pdf,
        ch_reports,
        multiqc_html
      )
    }
  emit:
    somatic_rds = ch_somatic_rds
    qc_pdf      = ch_qc_pdf
    tmb_csv     = ch_tmb_csv
    reports     = ch_reports
}

// --- Workflow: from_pretty_csv (downstream only) ---
workflow from_pretty_csv {
  main:
    ch_no_multiqc = Channel.fromPath("${projectDir}/assets/empty.txt").first()
    PARSE_PRETTY_CSV(ch_samplesheet_pretty)
    downstream_from_parse(PARSE_PRETTY_CSV.out, ch_no_multiqc)
}

// --- Workflow: from_fastq (upstream only: FASTQ -> pretty CSV) ---
workflow from_fastq {
  main:
    ch_ref_fasta = Channel.fromPath(params.ref_fasta).first()
    ch_dbsnp    = Channel.fromPath(params.dbsnp_vcf).first()
    ch_cosmic   = Channel.fromPath(params.cosmic_vcf).first()
    BUILD_BWA_INDEX(ch_ref_fasta)

    // FastQC on raw FASTQs (pre-trim)
    if (!params.skip_fastqc) {
      FASTQC_RAW(ch_fastq_rows.map { sid, r1, r2, pid, kit, proj, ar1, ar2 -> [ sid, r1, r2 ] })
    }

    TRIM_FASTQ(ch_fastq_rows)

    // FastQC on trimmed FASTQs (post-trim)
    if (!params.skip_fastqc) {
      FASTQC_TRIMMED(TRIM_FASTQ.out.map { sid, r1, r2, pid, kit, proj -> [ sid, r1, r2 ] })
    }

    ALIGN_BWA(TRIM_FASTQ.out, BUILD_BWA_INDEX.out.ref_index)
    ADD_READ_GROUPS(ALIGN_BWA.out)
    MARK_DUPLICATES(ADD_READ_GROUPS.out)
    BQSR(MARK_DUPLICATES.out.bam, ch_ref_fasta, ch_dbsnp)
    SAMTOOLS_STATS(BQSR.out)
    MUTECT(BQSR.out, ch_ref_fasta, ch_dbsnp, ch_cosmic)
    VARSCAN(BQSR.out, ch_ref_fasta)
    BCFTOOLS_STATS_MUTECT(MUTECT.out.map { sid, vcf, pid, kit, proj -> [ sid, vcf ] })
    BCFTOOLS_STATS_VARSCAN(VARSCAN.out.map { sid, snp, indel, pid, kit, proj -> [ sid, snp ] })
    // Join MUTECT (sample_id, vcf, patient_id, kit, project) with
    // VARSCAN (sample_id, snp_vcf, indel_vcf, patient_id, kit, project)
    // → merged tuple: (sample_id, mutect_vcf, varscan_snp_vcf, varscan_indel_vcf, patient_id, kit, project)
    ch_merge = MUTECT.out.join(VARSCAN.out).map { m, v -> [ m[0], m[1], v[1], v[2], m[2], m[3], m[4] ] }
    MERGE_CALLERS_ANNOVAR(ch_merge)
    COLLECT_PRETTY(MERGE_CALLERS_ANNOVAR.out)

    // Aggregate all QC with MultiQC
    ch_multiqc_html = Channel.fromPath("${projectDir}/assets/empty.txt").first()
    if (!params.skip_multiqc) {
      ch_qc = Channel.empty()
      if (!params.skip_fastqc) {
        ch_qc = ch_qc.mix(FASTQC_RAW.out.zip)
        ch_qc = ch_qc.mix(FASTQC_TRIMMED.out.zip)
      }
      ch_qc = ch_qc.mix(MARK_DUPLICATES.out.metrics)
      ch_qc = ch_qc.mix(SAMTOOLS_STATS.out.stats)
      ch_qc = ch_qc.mix(SAMTOOLS_STATS.out.flagstat)
      ch_qc = ch_qc.mix(SAMTOOLS_STATS.out.idxstats)
      ch_qc = ch_qc.mix(BCFTOOLS_STATS_MUTECT.out.stats)
      ch_qc = ch_qc.mix(BCFTOOLS_STATS_VARSCAN.out.stats)
      MULTIQC(ch_qc.collect())
      ch_multiqc_html = MULTIQC.out.report
    }

  emit:
    pretty_csvs  = COLLECT_PRETTY.out
    multiqc_html = ch_multiqc_html
}
