#!/usr/bin/env nextflow

/*
* Generate an experiment summary report containing statistics for multiple samples
*/
process multi_sample_report {
  publishDir "${params.outdir}/report/", mode: 'copy'

  input:
  path(csvs)
  path(template)
  path(qc_cascade_html)

  output:
  path('multisample_report.html')
  path('multisample_out.csv')
  path('multisample_summary_plots.html')
  path('multisample_qc_cascade.html')

  script:
  """
  create_multi_sample_report.py $template ${params.mixed_species} $qc_cascade_html
  """

  stub:
  """
  touch multisample_report.html
  touch multisample_out.csv
  touch multisample_summary_plots.html
  touch multisample_qc_cascade.html
  """
}
