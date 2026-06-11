#!/usr/bin/env nextflow

/*
* Generate QC cascade plot for single sample
*/
process qc_cascade_plot_single {
  tag "$sample_id"

  publishDir { "${params.outdir}/report/${sample_id}" }, mode: 'copy'

  input:
  tuple val(sample_id), path(metrics_csv)

  output:
  tuple val(sample_id), path("${sample_id}.qc_cascade.html"), emit: qc_cascade_plot

  script:
  """
  qc_cascade_plot.py \\
    --mode single \\
    --sample-id ${sample_id} \\
    --metrics-csv ${metrics_csv}
  """
}
