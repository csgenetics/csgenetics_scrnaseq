#!/usr/bin/env nextflow

/*
* Fused multimapper gene-assignment (WALL-CLOCK optimisation).
*
* Replaces four separate processes that ran as a per-sample serial chain, each paying a container
* start + S3 stage-in/out + queue wait on the cluster:
*   filter_for_multimappers_mismatch -> multimapper_transcript_assignment ->
*   multimapper_exon_assignment -> merge_transcript_exon_multimapper_bams
* The intermediate BAMs are now local to one task. The commands are byte-for-byte the same as the
* unfused chain (verified), so the output (${sample_id}.multimapped.annotated.bam) is unchanged.
*
* NB: the gawk (assign_multi_mappers.gawk) writes fixed filenames assigned_reads.sam_body /
* ambiguous_reads.sam_body and only (re)creates a file when it has something to write. In separate
* processes each got a fresh work dir; here we must rm them between the two gawk passes so a pass
* that produces no output does not leave a stale file from the previous pass.
*/
process multimapper_assignment {
  tag "$sample_id"

  input:
  tuple val(sample_id), path(feature_count_bam)
  path(gtf)
  path(multi_mapper_script)

  output:
  tuple val(sample_id), path("${sample_id}.multimapped.annotated.bam"), emit: high_conf_annotated_multimapped_bam

  script:
  """
  # 1. filter_for_multimappers_mismatch: keep multimappers (NH>1) with <=3 mismatches
  samtools view -@ ${task.cpus} -h -b -e '[NH]>1 && ([nM]==0 || [nM]==1 || [nM]==2 || [nM]==3)' -b ${feature_count_bam} > mm.featureCounts.bam

  # 2. multimapper_transcript_assignment: name-sort, gawk-split into assigned/ambiguous, rebuild bams
  samtools sort -n -@ ${task.cpus} -m 1G mm.featureCounts.bam | samtools view -@ ${task.cpus} | gawk -f ${multi_mapper_script}
  if [ -f assigned_reads.sam_body ]; then
    cat <(samtools view -H mm.featureCounts.bam) assigned_reads.sam_body | samtools view -@ ${task.cpus} -b -h > transcript.assigned.bam
  else
    samtools view -H -b mm.featureCounts.bam > transcript.assigned.bam
  fi
  if [ -f ambiguous_reads.sam_body ]; then
    cat <(samtools view -H mm.featureCounts.bam) ambiguous_reads.sam_body | samtools view -@ ${task.cpus} -h -b > transcript.unassigned.bam
  else
    samtools view -H -b mm.featureCounts.bam > transcript.unassigned.bam
  fi
  rm -f assigned_reads.sam_body ambiguous_reads.sam_body

  # 3. multimapper_exon_assignment: exon tie-break on the still-ambiguous reads, gawk-split again
  featureCounts -a ${gtf} -o exon.txt -R BAM transcript.unassigned.bam -T ${task.cpus} -t exon -g gene_id --fracOverlap 0.5 --extraAttributes gene_name -s 1 -M
  samtools sort -n -@ ${task.cpus} -m 1G transcript.unassigned.bam.featureCounts.bam | samtools view -@ ${task.cpus} | gawk -f ${multi_mapper_script}
  if [ -f assigned_reads.sam_body ]; then
    cat <(samtools view -H transcript.unassigned.bam.featureCounts.bam) assigned_reads.sam_body | samtools view -@ ${task.cpus} -b -h > exon.assigned.bam
  else
    samtools view -H -b transcript.unassigned.bam.featureCounts.bam > exon.assigned.bam
  fi

  # 4. merge the transcript- and exon-assigned bams
  samtools merge -@ ${task.cpus} -o ${sample_id}.multimapped.annotated.bam transcript.assigned.bam exon.assigned.bam
  """

  stub:
  """
  touch ${sample_id}.multimapped.annotated.bam
  """
}
