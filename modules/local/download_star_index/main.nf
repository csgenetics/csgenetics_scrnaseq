#!/usr/bin/env nextflow

/*
* Download the star index from the s3://csgx.public.readonly bucket
*/
process download_star_index {
  tag "Download STAR index"

  output:
  path("star"), emit: star_index

  script:
  """
  mkdir star
  aws s3 cp --no-sign-request ${params.star_index} ./star --recursive
  """
}
