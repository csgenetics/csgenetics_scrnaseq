#!/usr/bin/env nextflow

process CREATE_TEST_BAMS {
  container 'quay.io/biocontainers/samtools:1.17--hd87286a_1'

  input:
  path umr_sam
  path multimapper_sam

  output:
  path 'umr.bam', emit: umr
  path 'multimapper.bam', emit: multimapper

  script:
  """
  samtools view -b ${umr_sam} > umr.bam
  samtools view -b ${multimapper_sam} > multimapper.bam
  """
}
