#!/usr/bin/env nextflow

process run_rseqc {
  tag "$sample_id"
  
  publishDir "${params.outdir}/RSeQC/read_distribution", mode: 'copy', pattern: "*_rseqc_results.txt", saveAs: {"${sample_id}.${prefix}_RSeQC.txt"}

  input:
  tuple val(sample_id), path(bam), val(count)
  path(bed)
  path(empty_rseqc_template)
  val(prefix)

  output:
  tuple val(sample_id), path("*_rseqc_results.txt"), emit: rseqc_log

  script:
  """
  if [[ ${count} > 0 ]]
    then
      # rseqc_by_contig.py splits the BAM by reference contig and runs read_distribution.py on each
      # contig in parallel, then sums the per-feature counts (read_distribution is position-local, so
      # the output is byte-identical to running it on the whole BAM). Parallelises a single-threaded step.
      rseqc_by_contig.py -i ${bam} -r ${bed} -p ${task.cpus} > ${sample_id}_${prefix}_rseqc_results.txt
    else
      cat ${empty_rseqc_template} | envsubst > ${sample_id}_${prefix}_rseqc_results.txt
  fi
  """

  stub:
  """
  touch ${sample_id}_${prefix}_rseqc_results.txt
  """
}
