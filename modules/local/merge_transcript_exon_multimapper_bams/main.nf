#!/usr/bin/env nextflow

process merge_transcript_exon_multimapper_bams {
  tag "$sample_id"

  input:
  tuple val(sample_id), path(multimapper_transcript_bam), path(multimapper_exon_bam)

  output:
  tuple val(sample_id), path("${sample_id}.multimapped.annotated.bam"), emit: high_conf_annotated_multimapped_bam

  script:
  """
  samtools merge -o ${sample_id}.multimapped.annotated.bam $multimapper_transcript_bam $multimapper_exon_bam
  """
}
