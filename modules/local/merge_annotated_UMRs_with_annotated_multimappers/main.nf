#!/usr/bin/env nextflow

process merge_annotated_UMRs_with_annotated_multimappers {
  tag "$sample_id"

  publishDir "${params.outdir}/featureCounts", mode: 'copy', pattern: { "${sample_id}.annotated.bam" }

  input:
  tuple val(sample_id), path(umr_annotated_bam), path(multimapper_annotated_bam)

  output:
  tuple val(sample_id), path("${sample_id}.mapped.sorted.filtered.annotated.bam"), emit: high_conf_annotated_bam

  script:
  """
  samtools merge -o ${sample_id}.mapped.sorted.filtered.annotated.bam $umr_annotated_bam $multimapper_annotated_bam
  """
}
