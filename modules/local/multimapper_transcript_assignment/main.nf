#!/usr/bin/env nextflow

process multimapper_transcript_assignment{
  tag "$sample_id"

  input:
  tuple val(sample_id), path(multimapper_mismatch_filtered_bam)
  path(multi_mapper_script)

  output:
  tuple val(sample_id), path("${sample_id}.multimapped.transcript.assigned.bam"), emit: assigned_bam
  tuple val(sample_id), path("${sample_id}.multimapped.transcript.unassigned_ambiguity.no_xs_tag.bam"), emit: unassigned_bam

  script:
  """
  # Sort the bam file by query name in preparation for running through the gawk program for the first time.
  # Run the multimapped, annotated alignments through the gawk program that pulls
  # out alignments that can be associated directly as 'Assigned' and alignments that are ambiguous and should be passed onto
  # exon tie breaking.
  # See the script's header for more information on how it works.
  samtools sort -n $multimapper_mismatch_filtered_bam | samtools view | gawk -f $multi_mapper_script

  # multi_mapper_script produces assigned_reads.sam_body and ambiguous_reads.sam_body corresponding to the Assigned and still ambigous reads, respectively.
  # These files will only be produced if there were reads of the respective type identified.
  # If they exist, we need to convert both of these files back into valid bam format by adding a header and converting from sam->bam
  # If they don't exist then we should create a valid empty bam by the same name.

  if [ -f assigned_reads.sam_body ]; then
    # Cat with the headers of the featureCounts bam
    cat <(samtools view -H $multimapper_mismatch_filtered_bam) assigned_reads.sam_body | samtools view -b -h > ${sample_id}.multimapped.transcript.assigned.bam
  else
    samtools view -H -b $multimapper_mismatch_filtered_bam > ${sample_id}.multimapped.transcript.assigned.bam
  fi

  # If the ambigous_reads.sam_body exists then we convert this back to a valid bam.
  if [ -f ambiguous_reads.sam_body ]; then
    # Cat with the headers of the featureCounts bam
    # before rerunning through featureCounts for exon tie-breaking
    cat <(samtools view -H $multimapper_mismatch_filtered_bam) ambiguous_reads.sam_body | samtools view -h -b > ${sample_id}.multimapped.transcript.unassigned_ambiguity.no_xs_tag.bam
  else
    samtools view -H -b $multimapper_mismatch_filtered_bam > ${sample_id}.multimapped.transcript.unassigned_ambiguity.no_xs_tag.bam
  fi
  """

  stub:
  """
  touch ${sample_id}.multimapped.transcript.assigned.bam
  touch ${sample_id}.multimapped.transcript.unassigned_ambiguity.no_xs_tag.bam
  """
}
