#!/usr/bin/env nextflow

/*
* Categorize reads based on whether they are associated with cellular or non-cellular barcodes
* Also calculate counts in cells vs out of cells
*/
process categorize_reads {
  tag "$sample_id"

  input:
  tuple val(sample_id), path(star_bam), path(raw_count_matrix_h5ad), path(fastp_json)

  output:
  tuple val(sample_id), path("${sample_id}.read_categorization.csv"), emit: read_categories

  script:
  def barcode_length = params.barcode_pattern.count('C')
  def mixed_flag = params.mixed_species ? "--mixed_species" : ""
  """
  categorize_reads.py \
    --sample_id ${sample_id} \
    --star_bam ${star_bam} \
    --raw_count_matrix_h5ad ${raw_count_matrix_h5ad} \
    --fastp_json ${fastp_json} \
    --barcode_length ${barcode_length} \
    ${mixed_flag}
  """

  stub:
  """
  touch ${sample_id}.read_categorization.csv
  """
}
