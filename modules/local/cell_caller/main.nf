#!/usr/bin/env nextflow

/*
* Run cell caller - this determines a threshold number of counts to call a cell.
* The threshold is output to stdout and output as cell_caller_out
* If the h5ad matrix is empty, an empty .html will be output and checked for
* in the summary statistic script causing the Cell Caller plot to be hidden.
*/
process cell_caller {
  tag "$sample_id"

  // Publish the plots; the glob pattern is used to collect the mixed species and single species plots
  // Single species are named: {self.sample_name}_pdf_with_cutoff.html
  // Mixed species are named: {self.sample_name}_hsap_pdf_with_cutoff.html and {self.sample_name}_mmus_pdf_with_cutoff.html
  // pattern must be a string glob, not a closure (a closure stringifies and matches nothing on NF26).
  // The task is per-sample so a wildcard is unambiguous; matches both single (_pdf_with_cutoff.html)
  // and mixed (_hsap_/_mmus_pdf_with_cutoff.html) plots.
  publishDir "${params.outdir}/plots", pattern: "*_pdf_with_cutoff.html", mode: 'copy'

  input:
  tuple val(sample_id), path(count_matrix_h5ad), val(manual_threshold_str)

  output:
  tuple val(sample_id), stdout, emit: cell_caller_out
  tuple val(sample_id), path("${sample_id}_counts_pdf_with_threshold.html"), path("${sample_id}_barnyard_plot.html"), emit: cell_caller_plots
  // publishDir statement only works on files output in the output directive
  tuple val(sample_id), path("${sample_id}*_pdf_with_cutoff.html"), emit: cell_caller_html

  script:
  """
  cell_caller.py --sample_name ${sample_id} --minimum_count_threshold ${params.minimum_count_threshold} --count_matrix ${count_matrix_h5ad} --single_species ${!params.mixed_species} --manual_threshold_str $manual_threshold_str
  """

  stub:
  """
  touch ${sample_id}_counts_pdf_with_threshold.html
  touch ${sample_id}_barnyard_plot.html
  touch ${sample_id}_pdf_with_cutoff.html
  echo "100"
  """
}
