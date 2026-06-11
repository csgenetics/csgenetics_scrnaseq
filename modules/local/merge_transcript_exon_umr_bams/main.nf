#!/usr/bin/env nextflow

process merge_transcript_exon_umr_bams {
  tag "$sample_id"

  input:
  tuple val(sample_id), path(umr_transcript_bam), path(umr_exon_bam)

  output:
  tuple val(sample_id), path("${sample_id}.umr.annotated.bam"), emit: high_conf_annotated_umr_bam

  script:
  """
  samtools merge -o ${sample_id}.umr.annotated.bam $umr_transcript_bam $umr_exon_bam
  """

  stub:
  """
  touch ${sample_id}.umr.annotated.bam
  """
}
