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
  # alternative to umi_tools count
  # command: View the input dedup.sam file, grep for reads with gene assignment (which have the 'XT:' tag)
  # then keep the 1st (read_ID) and 18th fields (gene name, should be the 'XT:' tag itself) only
  # then split the string by '_' (-d '_'), keep fields 2,3,4 which correspond to 'io_sequence', 'gene_name' (4th field is so that we don't lose gene_names containing one '_')
  # use sed to replace the first '_' in each line, and any 'XT:Z:' strings with empty string with sed

  # output has 2 columns: io_sequence and gene_name for every deduplicated alignments with gene assignment
  samtools view ${f} | awk '/XT:/ {match(\$1, /_[A-Z]+_\$/); printf substr(\$0,RSTART+1,RLENGTH-2); match(\$0, /XT:Z:[A-Za-z0-9_]+/); print "\\t" substr(\$0,RSTART+5,RLENGTH-5)}' > ${sample_id}_bcGeneSummary.txt
  """

}
