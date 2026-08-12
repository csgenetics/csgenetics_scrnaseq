#!/usr/bin/env nextflow

process merge_annotated_UMRs_with_annotated_multimappers {
  tag "$sample_id"

  // pattern must be a string glob (a closure matches nothing on NF26) AND must match the actual output
  // filename below (${sample_id}.mapped.sorted.filtered.annotated.bam) -- the old glob matched neither.
  publishDir "${params.outdir}/featureCounts", mode: 'copy', pattern: "*.mapped.sorted.filtered.annotated.bam"

  input:
  tuple val(sample_id), path(umr_annotated_bam), path(multimapper_annotated_bam)

  output:
  tuple val(sample_id), path("${sample_id}.mapped.sorted.filtered.annotated.bam"), emit: high_conf_annotated_bam

  script:
  """
  # Both assignment branches may be name-ordered while retaining an upstream
  # SO:coordinate header. Canonically sort the merged customer-facing BAM so
  # its filename/header contract is true and it can be indexed directly.
  samtools merge -u -@ ${task.cpus} -o - \
    $umr_annotated_bam $multimapper_annotated_bam \
    | samtools sort -@ ${task.cpus} -m 1G \
        -o ${sample_id}.mapped.sorted.filtered.annotated.bam -
  """

  stub:
  """
  touch ${sample_id}.mapped.sorted.filtered.annotated.bam
  """
}
