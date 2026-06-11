#!/usr/bin/env nextflow

/*
* Convert GTF into features file
*/
process features_file {
  input:
  path(gtf)

  output:
  path("${gtf.baseName}_features_names.tsv"), emit: modified_gtf

  script:
  """
  features_names.py ${gtf} ${gtf.baseName}_features_names.tsv
  """

  stub:
  """
  touch ${gtf.baseName}_features_names.tsv
  """
}
