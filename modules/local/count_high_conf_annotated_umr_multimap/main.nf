#!/usr/bin/env nextflow

process count_high_conf_annotated_umr_multimap {
  tag "$sample_id"

  input:
  tuple val(sample_id), path(umr_multimapper_annotated_bam)

  output:
  tuple val(sample_id), env('alignment_count'), emit: aligned_count

  script:
  """
  alignment_count=\$(samtools view -c ${umr_multimapper_annotated_bam})
  """
}
