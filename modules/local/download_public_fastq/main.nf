#!/usr/bin/env nextflow

/* Download a fastq hosted on the s3://csgx.public.readonly bucket */
process download_public_fastq {
  tag "${sample_id} ${fastq_1.tokenize('/').last()}"

  input:
  tuple val(sample_id), val(fastq_1), val(fastq_2)

  output:
  tuple val(sample_id), path("${fastq_1.tokenize('/').last()}"), path("${fastq_2.tokenize('/').last()}"), emit: downloaded_fastqs

  script:
  def fastq_1_name = fastq_1.tokenize('/').last()
  def fastq_2_name = fastq_2.tokenize('/').last()
  """
  echo $fastq_1_name
  aws s3 cp --no-sign-request $fastq_1 $fastq_1_name
  aws s3 cp --no-sign-request $fastq_2 $fastq_2_name
  """
}
