#!/usr/bin/env nextflow

/*
* Filter count table - produce a count table that has been filtered for barcodes that pass the cell caller threshold
* If an empty file is input as the .h5ad (due to fail at count_matrix), an "${sample_id}.*.cell_only.count_matrix.empty.h5ad"
* file will be created and carried into summary_metrics.py. The tripartite matrix files will not be ouput in this case.
*/
process filter_count_matrix{
  tag "$sample_id"

  publishDir { "${params.outdir}/count_matrix/filtered_feature_bc_matrix/${sample_id}/" }, mode: 'copy', pattern: "*.filtered_feature_bc_matrix.h5ad"
  publishDir { "${params.outdir}/count_matrix/filtered_feature_bc_matrix/${sample_id}/" }, mode: 'copy', pattern: "matrix.mtx.gz"
  publishDir { "${params.outdir}/count_matrix/filtered_feature_bc_matrix/${sample_id}/" }, mode: 'copy', pattern: "barcodes.tsv.gz"
  publishDir { "${params.outdir}/count_matrix/filtered_feature_bc_matrix/${sample_id}/" }, mode: 'copy', pattern: "features.tsv.gz"
  publishDir { "${params.outdir}/count_matrix/raw_feature_bc_matrix/${sample_id}/" }, mode: 'copy', pattern: "*.raw_feature_bc_matrix.h5ad"

  input:
  tuple val(sample_id), val(count_threshold), path(h5ad_raw_count_matrix)

  output:
  // Output the 3 part barcode, features, matrix 
  tuple val(sample_id), path("barcodes.tsv.gz"), path("features.tsv.gz"), path("matrix.mtx.gz"), optional: true
  // Output the h5ad matrix
  tuple val(sample_id), path("${sample_id}*.filtered_feature_bc_matrix.*h5ad")
  // Output and emit the raw h5ad matrix for use in summary statistics
  tuple val(sample_id), path("${sample_id}*.raw_feature_bc_matrix.*h5ad"), emit: raw_count_matrix

  script:
  def mixed_args = params.mixed_species ? "TRUE" : "FALSE"
  """
  filter_count_matrix.py ${count_threshold} ${h5ad_raw_count_matrix} ${sample_id} $mixed_args
  """

  stub:
  // The raw h5ad is also an input (h5ad_raw_count_matrix) staged under the same
  // base name; Nextflow excludes staged inputs from output matching, so the stub
  // must emit a distinct, newly-created file that still matches the output glob.
  """
  touch barcodes.tsv.gz features.tsv.gz matrix.mtx.gz
  touch ${sample_id}.filtered_feature_bc_matrix.h5ad
  touch ${sample_id}.stub.raw_feature_bc_matrix.h5ad
  """
}
