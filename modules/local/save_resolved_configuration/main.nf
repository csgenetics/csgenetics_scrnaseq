#!/usr/bin/env nextflow

process save_resolved_configuration{
  tag "Save resolved configuration"
  publishDir "${params.outdir}/pipeline_info/", mode: 'copy'

  input:
  val(minimum_count_threshold)

  output:
  path("resolved_configuration.txt")

  script:
  // Params are immutable after workflow start. Replace the validated numeric
  // value in a copy so CLI/params-file strings cannot leak into provenance.
  resolved_params = new LinkedHashMap(params)
  resolved_params.minimum_count_threshold = minimum_count_threshold
  json_str = groovy.json.JsonOutput.toJson(resolved_params)
  json_indented = groovy.json.JsonOutput.prettyPrint(json_str)
  // Customer paths and CLI values may contain quotes, shell metacharacters or
  // Unicode. Carry the UTF-8 JSON through the generated shell only as base64.
  json_base64 = json_indented.getBytes('UTF-8').encodeBase64().toString()
  """
  printf '%s' '${json_base64}' | base64 -d > resolved_configuration.txt
  """

  stub:
  """
  touch resolved_configuration.txt
  """
}
