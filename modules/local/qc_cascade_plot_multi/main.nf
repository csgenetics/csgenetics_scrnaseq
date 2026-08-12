#!/usr/bin/env nextflow

/*
* Generate QC cascade plot for all samples (multi-sample)
*/
process qc_cascade_plot_multi {
  // The customer-facing plot embeds Plotly.js and works offline. Keep the
  // smaller fragment task-local for the consolidated report only.
  publishDir "${params.outdir}/report/", mode: 'copy', pattern: "multisample_qc_cascade.html"

  input:
  path(csvs)

  output:
  path("multisample_qc_cascade.fragment.html"), emit: qc_cascade_fragment
  path("multisample_qc_cascade.html"), emit: qc_cascade_report

  script:
  """
  qc_cascade_plot.py \\
    --mode multi \\
    --csv-files ${csvs}
  """

  stub:
  """
  touch multisample_qc_cascade.fragment.html
  touch multisample_qc_cascade.html
  """
}
