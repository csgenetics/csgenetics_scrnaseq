#!/usr/bin/env nextflow

/*
* Generate various count file precursors
*/
process io_count {
  tag "$sample_id"

  input:
  tuple val (sample_id), path(f)

  output:
  tuple val(sample_id), path('*_bcGeneSummary.txt'), emit: io_count_out

  script:
  """
  ## 1. Generate file for count matrix
  # For every deduplicated alignment carrying a gene assignment (the featureCounts 'XT:Z:' tag),
  # emit two columns: io_sequence (the barcode token at the end of the read name) and gene_name.
  #
  # io_count_extract is a static binary (tools/io_count_extract, on the pipeline bin/ PATH) that
  # replaces the previous awk one-liner while preserving complete featureCounts gene IDs and
  # normalising only terminal numeric version suffixes. The production samtools container ships
  # BusyBox awk, which is very slow (~409 s on a 1.4 GB BAM); the compiled, streaming replacement
  # performs the extraction in ~22 s. samtools still does the (fast, threaded) BAM decode.
  samtools view -@ ${task.cpus} ${f} | io_count_extract > ${sample_id}_bcGeneSummary.txt
  """

  stub:
  """
  touch ${sample_id}_bcGeneSummary.txt
  """

}
