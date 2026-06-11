#!/usr/bin/env nextflow

/*
* Generate a h5ad count matrix
* output raw_count_matrix (unfiltered for single-cells only)
* Will output ${sample_id}.count_matrix.empty.h5ad if the input file is blank
* Otherwise the file will be called ${sample_id}.count_matrix.h5ad
* The tripartite count table files will only be output if the input
* file is not empty.
*/
process count_matrix {
  tag "$sample_id"

  publishDir { "${params.outdir}/count_matrix/raw_feature_bc_matrix/${sample_id}/" }, mode: 'copy', pattern: "matrix.mtx.gz"
  publishDir { "${params.outdir}/count_matrix/raw_feature_bc_matrix/${sample_id}/" }, mode: 'copy', pattern: "barcodes.tsv.gz"
  publishDir { "${params.outdir}/count_matrix/raw_feature_bc_matrix/${sample_id}/" }, mode: 'copy', pattern: "features.tsv.gz"

  input:
  tuple val(sample_id), path(input_file)
  path(barcode_list)
  path(features_file)

  output:
  // Output the h5ad matrix
  tuple val(sample_id), path("${sample_id}.raw_feature_bc_matrix.*h5ad"), emit: h5ad
  // Output the 3 part barcode, features, matrix 
  tuple val(sample_id), path("barcodes.tsv.gz"), path("features.tsv.gz"), path("matrix.mtx.gz"), optional: true

  script:
  def mixed_args = params.mixed_species ? "--mixed_species True --hsap_mito_chr ${params.hsap_mitochondria_chromosome} --mmus_mito_chr ${params.mmus_mitochondria_chromosome} --hsap_gene_prefix ${params.hsap_gene_prefix} --mmus_gene_prefix ${params.mmus_gene_prefix}" : "--mixed_species False --mito_chr ${params.mitochondria_chromosome}"
  """
  count_matrix.py --barcode_list ${barcode_list} --count_table ${input_file} --gene_list ${features_file} --sample ${sample_id} $mixed_args
  """

  stub:
  """
  touch ${sample_id}.raw_feature_bc_matrix.h5ad
  touch barcodes.tsv.gz features.tsv.gz matrix.mtx.gz
  """
}
