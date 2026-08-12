#!/usr/bin/env nextflow

/*
* Unified QC process - replaces io_extract + io_extract_fastp + trim_extra_polya
* This process combines:
* 1. Barcode extraction and correction (io_extract functionality)
* 2. SSS trimming and polyX trimming (io_extract_fastp functionality)
* 3. Internal polyA trimming (trim_extra_polya functionality)
* 4. Q30 calculation for both R2 barcode and R1 output
* All in a single Rust binary for maximum performance
*/
process qc {
  tag "$sample_id"

  input:
  tuple val(sample_id), path(r1), path(r2)
  path(corrected_barcodelist)

  output:
  tuple val(sample_id), path("${sample_id}.qc.R1.fastq.gz"), emit: qc_out
  tuple val(sample_id), path("${sample_id}.qc.log"), emit: qc_log
  tuple val(sample_id), path("${sample_id}.R1.preQC.fastp.json"), path("${sample_id}.R2.preQC.fastp.json"), path("${sample_id}.R1.postQC.fastp.json"), emit: qc_multiqc

  script:
  """
  # Run unified qc Rust binary (available in PATH from container)
  qc \\
    ${r1} ${r2} ${corrected_barcodelist} \\
    --sample-id ${sample_id} \\
    --output-dir . \\
    --no-filtered

  # Check if output exists, create empty file if not
  if [ ! -f "${sample_id}.qc.R1.fastq.gz" ]; then
    touch ${sample_id}.qc.R1.fastq && gzip ${sample_id}.qc.R1.fastq
  fi
  """

  stub:
  def stubEmpty = (params.get('stub_empty_qc_samples') ?: []).contains(sample_id)
  def stubFastqCommand = stubEmpty ? "gzip -n </dev/null > ${sample_id}.qc.R1.fastq.gz" : "printf '@read1\\nACGT\\n+\\nIIII\\n' | gzip -n > ${sample_id}.qc.R1.fastq.gz"
  """
  # Tests can select QC-empty samples without changing production behavior.
  ${stubFastqCommand}
  touch ${sample_id}.qc.log
  touch ${sample_id}.R1.preQC.fastp.json ${sample_id}.R2.preQC.fastp.json ${sample_id}.R1.postQC.fastp.json
  """
}
