#!/usr/bin/env nextflow

// Filter alignments to allow for at most 3 mismatches to the reference.
process filter_for_UMRs_mismatch {
  tag "$sample_id"

  input:
  tuple val(sample_id), path(featurecount_bam)

  output:
  tuple val(sample_id), path("${sample_id}.UMRs.bam.featureCounts.bam"), emit: umr_mismatch_bam

  script:
  """
  samtools view -h -b -e '[NH]==1 && ([nM]==0 || [nM]==1 || [nM]==2 || [nM]==3)' -b ${featurecount_bam} > ${sample_id}.UMRs.bam.featureCounts.bam
  """
}
