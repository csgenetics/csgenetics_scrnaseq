#!/usr/bin/env nextflow

process save_resolved_configuration{
  tag "Save resolved configuration"
  publishDir "${params.outdir}/pipeline_info/", mode: 'copy'

  output:
  path("resolved_configuration.txt")

  script:
  json_str = groovy.json.JsonOutput.toJson(params)
  json_indented = groovy.json.JsonOutput.prettyPrint(json_str)
  // NOTE: single quotes are critical here;
  """
  echo '${json_indented}' > resolved_configuration.txt
  """
}
