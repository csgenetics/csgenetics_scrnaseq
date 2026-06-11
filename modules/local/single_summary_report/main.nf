#!/usr/bin/env nextflow

/*
* Generate a per sample html report
*/
process single_summary_report {
  tag "$sample_id"

  publishDir { "${params.outdir}/report/${sample_id}" }, mode: 'copy'

  input:
  tuple val(sample_id), path(metrics_csv), path(pdf_plot_html), path(barnyard_plot_html), path(qc_cascade_html)
  path(html_template)

  output:
  tuple val(sample_id), path("${sample_id}_report.html")
  path("${sample_id}.metrics.csv"), emit: single_sample_metric_out

  script:
  """
  create_single_sample_report.py $sample_id $pdf_plot_html $barnyard_plot_html $qc_cascade_html $metrics_csv $html_template ${params.mixed_species}
  """

  stub:
  """
  touch ${sample_id}_report.html
  touch ${sample_id}.metrics.csv
  """
}
