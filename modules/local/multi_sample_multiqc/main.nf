#!/usr/bin/env nextflow

/*
* Run multiqc to generate multisample report and data
*/
process multi_sample_multiqc {
  publishDir "${params.outdir}/multiqc/", mode: 'copy', pattern: "multisample_multiqc_data"
  publishDir "${params.outdir}/multiqc/", mode: 'copy', pattern: "multisample_multiqc.html"
  publishDir "${params.outdir}/multiqc/", mode: 'copy', pattern: "multisample_multiqc_data/multiqc_data.json", saveAs: {"multisample.multiqc.data.json"}

  input:
  path(multiqc_in_files)

  output:
  path "multisample_multiqc.html"
  path "multisample_multiqc_data"
  path("multisample_multiqc_data/multiqc_data.json")

  script:
  """
  multiqc . \
    -f \
    --title "multisample multiqc" \
    --filename "multisample_multiqc.html" \
    -m unified_qc \
    -m rseqc
  """
}
