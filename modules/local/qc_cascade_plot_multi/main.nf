#!/usr/bin/env nextflow

/*
* Generate QC cascade plot for all samples (multi-sample)
*/
process qc_cascade_plot_multi {
  publishDir "${params.outdir}/report/", mode: 'copy'

  input:
  path(csvs)

  output:
  path("multisample_qc_cascade.html"), emit: qc_cascade_plot

  script:
  """
  qc_cascade_plot.py \\
    --mode multi \\
    --csv-files ${csvs}
  """
}
