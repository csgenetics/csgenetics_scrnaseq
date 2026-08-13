#!/usr/bin/env nextflow

/* Download the gtf file from the s3://csgx.public.readonly bucket */
process download_gtf {
  tag "Download GTF"

  output:
  path("*.gtf")

  script:
  """
  aws s3 cp --no-sign-request ${params.gtf} .
  """

  stub:
  """
  touch genes.gtf
  """
}
