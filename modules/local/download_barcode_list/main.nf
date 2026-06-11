#!/usr/bin/env nextflow

/* Download the csgx hosted barcode csv
* from the s3://csgx.public.readonly bucket.
*/
process download_barcode_list {
  tag "Get barcode list"

  output:
  path("*.csv")

  script:
  """
  aws s3 cp --no-sign-request ${params.barcode_list_path} .
  """

  stub:
  """
  touch barcode_list.csv
  """
}
