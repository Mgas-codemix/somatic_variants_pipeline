// Build BWA index from reference FASTA. Used by from_fastq so alignment does not require a pre-built index.

process BUILD_BWA_INDEX {
  tag "ref"
  publishDir "${params.outdir}/reference", mode: params.publish_mode, pattern: "*", when: false  // optional: publish index for reuse
  input:
    path ref_fasta
  output:
    tuple path(ref_fasta), path("*.amb"), path("*.ann"), path("*.bwt"), path("*.pac"), path("*.sa"), emit: ref_index
  script:
    def ref = ref_fasta.name
    def refPath = ref_fasta.toString()
    def algo = new File(refPath).length() < 200_000_000 ? "is" : "bwtsw"
    """
    ${params.bwa} index -a ${algo} ${ref}
    """
}
