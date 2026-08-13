#!/usr/bin/env nextflow

/* Can also be used to download the user specified cell caller threshold csv, if specified */
process download_input_csv {
  tag "Download input csvs"

  output:
  path("*.csv")

  script:
  """
  aws s3 cp --no-sign-request ${params.input_csv} .
  """

  stub:
  """
  touch input.csv
  """
}
