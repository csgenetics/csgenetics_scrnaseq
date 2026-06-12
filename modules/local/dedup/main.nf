#!/usr/bin/env nextflow

/*
* Deduplicate reads (currently based on the combination of IO and SSS start site)
*/
process dedup{
  tag "$sample_id"

  publishDir "${params.outdir}/deduplication", mode: 'copy'

  input:
  tuple val (sample_id), path(bam), path(bai), val(alignment_count)

  output:
  tuple val(sample_id), path("${sample_id}.dedup.log"), emit: io_dedup_log
  tuple val(sample_id), path("${sample_id}.dedup.bam"), emit: io_dedup_sam
  
  script:
  """
  # The effect of this will be to deduplicate based on the combination of cell
  # barcode and SSS start site.
  if [[ $alignment_count > 0 ]]
    then
      # umi_tools dedup is single-threaded and position-local. Split by reference contig, dedup
      # each contig in parallel, and merge (count-identical to whole-BAM dedup; see the script).
      dedup_by_contig.sh ${bam} ${sample_id} ${task.cpus}
    else
    
    # Then there were no alignments and we should output a dummy .dedup.log
    # and the original bam renamed
      echo "INFO Reads: Input Reads: 0\nINFO Number of reads out: 0\n" > ${sample_id}.dedup.log
      cp $bam ${sample_id}.dedup.bam
  fi
  """

  stub:
  """
  touch ${sample_id}.dedup.log
  touch ${sample_id}.dedup.bam
  """
}
