#!/usr/bin/env nextflow

/*
* Generate the summary statistics required for HTML report
*/
process summary_statistics {
  tag "$sample_id"

  publishDir { "${params.outdir}/report/${sample_id}" }, mode: 'copy', pattern: "*.csv"

  input:
  tuple val(sample_id), val(minimum_count_threshold), path(raw_h5ad), path(antisense), path(dedup), path("${sample_id}.multiqc.data.json"), path("${sample_id}_raw_rseqc_results.txt"), path("${sample_id}_annotated_rseqc_results.txt"), path(read_categorization_csv), path(qc_log)
  output:
  tuple val(sample_id), path("${sample_id}.metrics.csv"), emit: metrics_csv

  script:
  def mixed_flag = params.mixed_species ? "--mixed-species" : ""
  """
  summary_statistics.py \
    --sample-id ${sample_id} \
    --h5ad ${raw_h5ad} \
    --multiqc-json ${sample_id}.multiqc.data.json \
    --antisense ${antisense} \
    --dedup-log ${dedup} \
    --raw-rseqc ${sample_id}_raw_rseqc_results.txt \
    --annotated-rseqc ${sample_id}_annotated_rseqc_results.txt \
    --read-categorization ${read_categorization_csv} \
    --qc-log ${qc_log} \
    ${mixed_flag}
  """

  stub:
  """
  touch ${sample_id}.metrics.csv
  """
}
