#!/usr/bin/env nextflow

process run_rseqc {
  tag "$sample_id"
  
  publishDir "${params.outdir}/RSeQC/read_distribution", mode: 'copy', pattern: "*_rseqc_results.txt", saveAs: {"${sample_id}.${prefix}_RSeQC.txt"}

  input:
  tuple val(sample_id), path(bam), val(alignment_count)
  path(bed)
  path(empty_rseqc_template)
  val(prefix)

  output:
  tuple val(sample_id), path("*_rseqc_results.txt"), emit: rseqc_log

  script:
  """
  if [[ ! "${alignment_count}" =~ ^[0-9]+\$ ]]
    then
      echo "alignment_count must be a non-negative integer; got '${alignment_count}'" >&2
      exit 1
  elif [[ ! "${alignment_count}" =~ ^0+\$ ]]
    then
      # NB: a contig-split parallel variant (bin/rseqc_by_contig.py) was tried but is SLOWER on real
      # large BAMs -- the single-threaded pysam split (one full read pass) costs more than the
      # per-contig parallelism saves (raw_rseqc 6.5 -> 25.6 min/sample). Reverted to the direct call.
      # A future speedup would need a samtools-based scatter (the rseqc container has no samtools).
      read_distribution.py -i ${bam} -r ${bed} > ${sample_id}_${prefix}_rseqc_results.txt
    else
      cat ${empty_rseqc_template} | envsubst > ${sample_id}_${prefix}_rseqc_results.txt
  fi
  """

  stub:
  """
  if [[ ! "${alignment_count}" =~ ^[0-9]+\$ ]]
    then
      echo "alignment_count must be a non-negative integer; got '${alignment_count}'" >&2
      exit 1
  elif [[ ! "${alignment_count}" =~ ^0+\$ ]]
    then
      printf 'Total Tags                    1\n' > ${sample_id}_${prefix}_rseqc_results.txt
    else
      printf 'Total Tags                    0\n' > ${sample_id}_${prefix}_rseqc_results.txt
  fi
  """
}
