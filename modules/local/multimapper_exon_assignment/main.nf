#!/usr/bin/env nextflow

process multimapper_exon_assignment{
  tag "$sample_id"

  input:
  tuple val(sample_id), path(multimapper_unassigned_bam)
  path(gtf)
  path(multi_mapper_script)

  output:
  tuple val(sample_id), path("${sample_id}.multimapped.exon.assigned.bam"), emit: assigned_bam

  script:
  """
  # Run featureCounts using the exon feature to do tie-breaking and then run back through the gawk script to pull out those
  # reads that have a single Assigned alignment.
  featureCounts -a $gtf -o ${sample_id}.multimapped.exon.assigned.txt -R BAM $multimapper_unassigned_bam -T ${task.cpus} -t exon -g gene_id --fracOverlap 0.5 --extraAttributes gene_name -s 1 -M
  # Threaded name sort (-@); the gawk is order-independent + deterministic (see assign_multi_mappers.gawk).
  samtools sort -n -@ ${task.cpus} -m 1G ${sample_id}.multimapped.transcript.unassigned_ambiguity.no_xs_tag.bam.featureCounts.bam | samtools view -@ ${task.cpus} | gawk -f $multi_mapper_script

  # If the assigned_reads.sam_body file exists then we were successfuly able to pull out further assigned reads
  if [ -f assigned_reads.sam_body ]; then
    # Cat with the headers of the featureCounts bam
    cat <(samtools view -H ${sample_id}.multimapped.transcript.unassigned_ambiguity.no_xs_tag.bam.featureCounts.bam) assigned_reads.sam_body | samtools view -@ ${task.cpus} -b -h > ${sample_id}.multimapped.exon.assigned.bam;
  else
    # There were no further reads successfuly annotated
    # Create a valid empty bam
    samtools view -H -b ${sample_id}.multimapped.transcript.unassigned_ambiguity.no_xs_tag.bam.featureCounts.bam > ${sample_id}.multimapped.exon.assigned.bam
  fi
  """

  stub:
  """
  touch ${sample_id}.multimapped.exon.assigned.bam
  """
}
