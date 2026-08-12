#!/usr/bin/env nextflow

/*
* Generate ONE consolidated, self-contained experiment report.
*
* Replaces the previous per-sample reports (${sample_id}_report.html) and the
* multi-sample report family (multisample_report.html, multisample_summary_plots.html).
* The per-sample ${sample_id}.metrics.csv (published by summary_statistics) and
* multisample_out.csv (re-emitted here, unchanged) are preserved.
*
* All per-sample metrics csvs and plot fragments are staged FLAT; the script
* discovers samples from the *.metrics.csv files. Vendored assets (Plotly.js,
* Bootstrap, fonts) are embedded inline so the report is offline-safe.
*/
process consolidated_report {
  publishDir "${params.outdir}/report/", mode: 'copy'

  input:
  path(metrics_csvs)
  path(cell_caller_plots)
  path(qc_cascade_plots)
  path(multi_qc_cascade_html)
  path(template)
  path(vendor_dir)
  val(provenance_base64)

  output:
  path('consolidated_report.html')
  path('multisample_out.csv')

  script:
  // Groovy base64-encodes run metadata before it reaches this process. The
  // resulting alphabet contains no quotes or shell metacharacters, so customer
  // paths and run names remain inert data across Nextflow's generated wrapper.
  """
  create_consolidated_report.py \\
    ${template} \\
    ${params.mixed_species} \\
    ${vendor_dir} \\
    ${multi_qc_cascade_html} \\
    'base64:${provenance_base64}'
  """

  stub:
  """
  touch consolidated_report.html
  touch multisample_out.csv
  """
}
