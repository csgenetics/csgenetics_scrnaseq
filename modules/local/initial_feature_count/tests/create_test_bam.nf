#!/usr/bin/env nextflow

process CREATE_TEST_BAM {
  container 'quay.io/biocontainers/samtools:1.17--hd87286a_1'

  input:
  path sam

  output:
  path 'aligned.bam', emit: bam

  script:
  """
  samtools view -b ${sam} > aligned.bam
  """
}
