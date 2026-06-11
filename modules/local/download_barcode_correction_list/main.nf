#!/usr/bin/env nextflow

/* Download the csgx hosted barcode correction list tsv
* from the s3://csgx.public.readonly bucket.
*/
process download_barcode_correction_list {
  tag "Get barcode correction list"

  output:
  path("*.tsv")

  script:
  """
  aws s3 cp --no-sign-request ${params.barcode_correction_list_path} .
  """
}
