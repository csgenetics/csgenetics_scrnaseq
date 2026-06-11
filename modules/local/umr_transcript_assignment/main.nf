#!/usr/bin/env nextflow

// Collect those UMRs that have unambiguous transcript annotations based on transcript feature (i.e. no exon tie-breaking required)
process umr_transcript_assignment {
  tag "$sample_id"

  input:
  tuple val(sample_id), path(umr_mismatch_bam)

  output:
  tuple val(sample_id), path("${sample_id}.UMRs.transcript.assigned.bam"), emit: umr_transcript_assigned_bam

  script:
  """
  # Filter for the reads that were 'Assigned' a transcript target
  samtools view -h -b -e '[XN]==1 && [XT] && [XS]=="Assigned"' -b ${umr_mismatch_bam} > ${sample_id}.UMRs.transcript.assigned.bam
  """

  stub:
  """
  touch ${sample_id}.UMRs.transcript.assigned.bam
  """
}
