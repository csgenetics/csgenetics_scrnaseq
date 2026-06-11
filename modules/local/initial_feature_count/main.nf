#!/usr/bin/env nextflow

process initial_feature_count {
  tag "$sample_id"

  publishDir "${params.outdir}/featureCounts", mode: 'copy'

  input:
  tuple val(sample_id), path(bam), val(aligned_count)
  path(gtf)

  output:
  tuple val(sample_id), path("${sample_id}_Aligned.sortedByCoord.out.bam.featureCounts.bam"), emit: feature_count_bam

  script:
  """
  if [[ $aligned_count > 0 ]] # If the bam is not empty
  then
    samtools sort $bam -o ${sample_id}_Aligned.sortedByCoord.out.bam
    # Start by running feature counts on the star output
    # including strandedness and annotation of multimappers
    featureCounts -a $gtf -o ${sample_id}.star.featureCounts.gene.txt -R BAM ${sample_id}_Aligned.sortedByCoord.out.bam -T ${task.cpus} -t transcript -g gene_id --fracOverlap 0.5 --extraAttributes gene_name -s 1 -M
  else
    # Simply rename the input bam so that it can be collected
    cp $bam ${sample_id}_Aligned.sortedByCoord.out.bam.featureCounts.bam
  fi
  """

  stub:
  """
  touch ${sample_id}_Aligned.sortedByCoord.out.bam.featureCounts.bam
  """
}
