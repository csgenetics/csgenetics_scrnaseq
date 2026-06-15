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
  # Name-hash-parallel: split the multimapper BAM by read-name hash, run the filter+transcript+exon+merge
  # sub-pipeline on each chunk in parallel, and merge. Parallelises the two name sorts (the dominant cost)
  # bounded by the largest chunk. Read-set identical to the single-pass (verified). See the script header.
  multimapper_assignment.sh ${feature_count_bam} ${sample_id} ${gtf} ${multi_mapper_script} ${task.cpus}
  """

  stub:
  """
  touch ${sample_id}.multimapped.annotated.bam
  """
}
