#!/usr/bin/env nextflow

process filter_for_multimappers_mismatch {
  tag "$sample_id"

  input:
  tuple val(sample_id), path(feature_count_bam)

  output:
  tuple val(sample_id), path("${sample_id}.multimapped.bam.featureCounts.bam"), emit: multimap_mismatch_bam

  script:
  """
  # -@ threads only the BGZF (de)compression; the per-read filter and output order are unchanged (byte-identical).
  samtools view -@ ${task.cpus} -h -b -e '[NH]>1 && ([nM]==0 || [nM]==1 || [nM]==2 || [nM]==3)' -b ${feature_count_bam} > ${sample_id}.multimapped.bam.featureCounts.bam
  """

  stub:
  """
  touch ${sample_id}.multimapped.bam.featureCounts.bam
  """
}
