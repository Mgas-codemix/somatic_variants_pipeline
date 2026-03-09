// Optional Cutadapt + SeqTK trimfq. Input: tuple (sample_id, fastq_r1, fastq_r2, patient_id, kit, project, adapter_r1, adapter_r2)
// Output: trimmed R1 and R2 (or single R1 for single-end)

process TRIM_FASTQ {
  tag "${sample_id}"
  input:
    tuple val(sample_id), path(fastq_r1), path(fastq_r2), val(patient_id), val(kit), val(project), val(adapter_r1), val(adapter_r2)
  output:
    tuple val(sample_id), path("*.trimmed.R1.fastq.gz"), path("*.trimmed.R2.fastq.gz"), val(patient_id), val(kit), val(project), emit: trimmed
  script:
    def r2 = fastq_r2.name != "null" && fastq_r2.name != "empty.txt"
    def ad1 = adapter_r1 ?: params.adapter_r1
    def ad2 = adapter_r2 ?: params.adapter_r2
    def skip_cut = params.skip_cutadapt ? "true" : "false"
    if (r2) {
      """
      if [ "${skip_cut}" != "true" ]; then
        cutadapt -a ${ad1} -A ${ad2} -m 20 -o ${sample_id}.clipped.R1.fastq.gz -p ${sample_id}.clipped.R2.fastq.gz ${fastq_r1} ${fastq_r2} || true
        IN1=${sample_id}.clipped.R1.fastq.gz; IN2=${sample_id}.clipped.R2.fastq.gz;
      else
        IN1=${fastq_r1}; IN2=${fastq_r2};
      fi
      seqtk trimfq -b 5 \$IN1 | gzip -c > ${sample_id}.trimmed.R1.fastq.gz
      seqtk trimfq -b 5 \$IN2 | gzip -c > ${sample_id}.trimmed.R2.fastq.gz
      """
    } else {
      """
      if [ "${skip_cut}" != "true" ]; then
        cutadapt -a ${ad1} -m 20 -o ${sample_id}.clipped.R1.fastq.gz ${fastq_r1} || true
        IN1=${sample_id}.clipped.R1.fastq.gz
      else
        IN1=${fastq_r1}
      fi
      seqtk trimfq -b 5 \$IN1 | gzip -c > ${sample_id}.trimmed.R1.fastq.gz
      touch ${sample_id}.trimmed.R2.fastq.gz
      """
    }
}
