#!/usr/bin/env nextflow

/*
* Sort and index bam file using samtools
*/
process sort_index_bam {
  tag "$sample_id"

  input:
  tuple val (sample_id), path(bam)

  output:
  tuple val(sample_id), path("${sample_id}.antisense.txt"), emit: antisense_out
  tuple val (sample_id), path("${sample_id}_sorted.bam"), path("${sample_id}_sorted.bam.bai"), emit: sort_index_bam_out

  script:
  """
  samtools view -c -f 16 ${bam} > ${sample_id}.antisense.txt
  samtools sort ${bam} -O BAM -o ${sample_id}_sorted.bam
  samtools index ${sample_id}_sorted.bam
  """

  stub:
  """
  touch ${sample_id}.antisense.txt
  touch ${sample_id}_sorted.bam
  touch ${sample_id}_sorted.bam.bai
  """
}
