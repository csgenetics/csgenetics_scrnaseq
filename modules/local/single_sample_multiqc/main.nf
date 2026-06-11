#!/usr/bin/env nextflow

// TODO Pick up the multiqc_data.json directly and rename it in summary_statistics
// NOTE I am going to see if we can pick this up without a double star glob.
process single_sample_multiqc {
  publishDir { "${params.outdir}/multiqc/single_sample_multiqc/${sample_id}" }, mode: 'copy', pattern: "*_data"
  publishDir { "${params.outdir}/multiqc/single_sample_multiqc/${sample_id}" }, mode: 'copy', pattern: "*_multiqc.html"
  publishDir { "${params.outdir}/multiqc/single_sample_multiqc/${sample_id}" }, mode: 'copy', pattern: "**/multiqc_data.json", saveAs: {"${sample_id}.multiqc.data.json"}

  input:
  tuple val(sample_id), path(multiqc_in_files)

  output:
  path "*_multiqc.html"
  path "*_data"
  tuple val(sample_id), path("**/multiqc_data.json"), emit: multiqc_json

  script:
  """
  multiqc . \
    -f \
    --title "${sample_id} multiqc" \
    --filename "${sample_id}_multiqc.html" \
    -m unified_qc \
    -m rseqc
  """

  stub:
  """
  touch ${sample_id}_multiqc.html
  mkdir -p ${sample_id}_multiqc_data
  touch ${sample_id}_multiqc_data/multiqc_data.json
  """
}
