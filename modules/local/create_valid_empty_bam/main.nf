#!/usr/bin/env nextflow

/*
* Create an empty bam file that has a single line header
* so that the downstream samtools filtering doesn't break.
* Otherwise STAR and featureCount produce an empty bam that is an empty file
* that breaks samtools view: [main_samview] fail to read the header from "247intergenic_Aligned.sortedByCoord.out.bam".
*/
process create_valid_empty_bam{
  tag "$sample_id"

  publishDir "${params.outdir}/STAR", mode: 'copy'

  input:
  tuple val(sample_id), val(prefix)

  output:
  tuple val(sample_id), path("${sample_id}${prefix}.bam"), emit: out_bam

  script:
  """
  echo "@HD	VN:1.4	SO:coordinate" | samtools view -h -b > ${sample_id}${prefix}.bam
  """

}
