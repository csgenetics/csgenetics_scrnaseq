#!/usr/bin/env nextflow

process umr_exon_assignment {
  tag "$sample_id"

  input:
  tuple val(sample_id), path(umr_mismatch_filtered_bam)
  path(gtf)

  output:
  tuple val(sample_id), path("${sample_id}.UMRs.exon.assigned.bam"), emit: umr_exon_assigned_bam

  script:
  """
  # Filter out the reads that were that were classified as Unassigned_Ambiguity and run them through
  # featureCounts using the exon tag to see if the ambiguity can be cleared up based on exon mapping.
  samtools view -h -b -e '[XS]=="Unassigned_Ambiguity"' -b $umr_mismatch_filtered_bam > ${sample_id}.UMRs.transcript.unassigned_ambiguity.bam

  # We have to remove the XS tag from the *.transcript.unassigned_ambiguity.bam because featureCounts adds
  # an additional tag, that prevents proper fitltering with samtools
  samtools view -h ${sample_id}.UMRs.transcript.unassigned_ambiguity.bam | sed 's/\\tXS\\:Z\\:[^\\t]*//' | samtools view -h -b > ${sample_id}.UMRs.transcript.unassigned_ambiguity.no_xs_tag.bam

  # Do exon tie-breaking and filter to those that are assigned.
  featureCounts -a $gtf -o ${sample_id}.UMRs.exon.assigned.txt -R BAM ${sample_id}.UMRs.transcript.unassigned_ambiguity.no_xs_tag.bam -T ${task.cpus} -t exon -g gene_id --fracOverlap 0.5 --extraAttributes gene_name -s 1
  samtools view -h -b -e '[XN]==1 && [XT] && [XS]=="Assigned"' -b ${sample_id}.UMRs.transcript.unassigned_ambiguity.no_xs_tag.bam.featureCounts.bam > ${sample_id}.UMRs.exon.assigned.bam
  """
}
