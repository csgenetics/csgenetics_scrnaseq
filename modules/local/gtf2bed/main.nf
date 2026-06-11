#!/usr/bin/env nextflow

/*
Process to convert the input GTF to a gene model bed file for rseqc read distribution
*/
process gtf2bed {
  input:
  path(gtf)

  output:
  path("gene_model.bed"), emit: bed

  script:
  """
  gtf2bed ${gtf} > gene_model.bed
  """

  stub:
  """
  touch gene_model.bed
  """
}
