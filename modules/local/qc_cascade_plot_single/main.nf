#!/usr/bin/env nextflow

/*
* Generate QC cascade plot for single sample
*/
process qc_cascade_plot_single {
  tag "$sample_id"

  // Only publish the self-contained output. The Plotly-free fragment is an
  // internal input to consolidated_report and would render blank by itself.
  publishDir { "${params.outdir}/report/${sample_id}" }, mode: 'copy', pattern: "*.qc_cascade.standalone.html", saveAs: { "${sample_id}.qc_cascade.html" }

  input:
  tuple val(sample_id), path(metrics_csv)

  output:
  tuple val(sample_id), path("${sample_id}.qc_cascade.html"), emit: qc_cascade_fragment
  tuple val(sample_id), path("${sample_id}.qc_cascade.standalone.html"), emit: qc_cascade_report

  script:
  """
  qc_cascade_plot.py \\
    --mode single \\
    --sample-id ${sample_id} \\
    --metrics-csv ${metrics_csv}
  """

  stub:
  """
  touch ${sample_id}.qc_cascade.html
  touch ${sample_id}.qc_cascade.standalone.html
  """
}
