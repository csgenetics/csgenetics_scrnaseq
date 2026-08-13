#!/usr/bin/env nextflow

process CREATE_EMPTY_BAM {
  container 'quay.io/biocontainers/samtools:1.17--hd87286a_1'

  output:
  path 'empty.bam', emit: bam

  script:
  """
  printf '@HD\tVN:1.4\tSO:coordinate\n' | samtools view -h -b > empty.bam
  """
}
