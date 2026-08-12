#!/usr/bin/env nextflow

/*
* Align reads with STAR
* We do the sorting separately in samtools to avoid the sorting RAM error
* thrown by STAR
*/
process star {
  tag "$sample_id"

  publishDir "${params.outdir}/STAR", mode: 'copy'

  input:
  tuple val(sample_id), path(r1)
  path(index)

  output:
  tuple val(sample_id), path("${sample_id}_Aligned.out.bam"), env('uniquely_mapped_reads'), env('aligned_reads'), emit: out_bam

  script:
  """
      # runThreadN is PINNED at 8 (not tied to task.cpus) on purpose. With --outFilterMultimapNmax
      # 1000 STAR's per-read multimapper output ORDER is thread-count-sensitive, so the thread count
      # is part of the OUTPUT CONTRACT: changing it shifts which multimappers get assigned and hence
      # per-barcode total counts. The original pipeline used `--runThreadN 8`; we keep 8 so alignment
      # output is identical to the original (a prior change to \${task.cpus}=16 silently altered the
      # counts). The cpu RESERVATION (conf/base.config) is right-sized independently for throughput.
      # Do NOT couple this back to \${task.cpus} and do NOT change the value - it changes results.
      STAR --runThreadN 8 \
        --genomeDir ${index} \
        --readFilesIn ${r1} \
        --outFileNamePrefix ${sample_id}_ \
        --outReadsUnmapped Fastx \
        --outSAMtype BAM Unsorted \
        --readFilesCommand zcat \
        --outSAMattributes Standard \
        --outFilterMultimapNmax 1000

      # STAR reports unique and accepted multi-locus mappings separately. Both
      # represent real aligned reads for downstream feature assignment. The
      # parser validates required fields, integer bounds, duplicates, overflow,
      # and consistency with the number of input reads before routing the sample.
      star_counts=\$(star_alignment_counts.py ${sample_id}_Log.final.out)
      IFS=\$'\t' read -r uniquely_mapped_reads aligned_reads <<< "\${star_counts}"

  """

  stub:
  def stubMultimapperOnly = (params.get('stub_multimapper_only_samples') ?: []).contains(sample_id)
  def stubUnaligned = (params.get('stub_unaligned_samples') ?: []).contains(sample_id)
  def stubUniqueCount = (stubMultimapperOnly || stubUnaligned) ? 0 : 1
  def stubAlignedCount = stubUnaligned ? 0 : 1
  """
  touch ${sample_id}_Aligned.out.bam
  uniquely_mapped_reads=${stubUniqueCount}
  aligned_reads=${stubAlignedCount}
  """

}
