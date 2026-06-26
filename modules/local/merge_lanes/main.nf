#!/usr/bin/env nextflow

/*
* Merge Lanes
*/
process merge_lanes {
  tag "$sample_id $read_num"

  input:
  tuple val(sample_id), val(read_num), path(fastqs)

  output:
  tuple val(sample_id), val(read_num), path("${sample_id}.merged.${read_num}.fastq.gz"), emit: merge_lanes_out

  script:
  """
  # cat ${fastqs} > ${sample_id}.merged.${read_num}.fastq.gz

  # Globs are ordered so lane merging will happen in same order for R1 and R2
  cat *.f*q.gz > ${sample_id}.merged.${read_num}.fastq.gz
  """

  stub:
  """
  touch ${sample_id}.merged.${read_num}.fastq.gz
  """
}
