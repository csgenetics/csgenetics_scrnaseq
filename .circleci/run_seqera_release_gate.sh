#!/usr/bin/env bash

# Launch the exact CircleCI revision on Seqera Platform, wait for a terminal
# result, and remove the temporary action. This script deliberately owns the
# complete lifecycle in one shell so EXIT/TERM cleanup can cancel an incomplete
# workflow and delete its action.

set -Eeuo pipefail

readonly seqera_api_base='https://api.cloud.seqera.io'
readonly pipeline_url='https://github.com/csgenetics/csgenetics_scrnaseq'
readonly nextflow_version='26.04.1'

action_id=''
workflow_id=''
workflow_terminal=0
launch_attempted=0
scratch_dir=''
auth_header_file=''

error() {
  printf 'ERROR: %s\n' "$1" >&2
}

die() {
  error "$1"
  exit 1
}

require_environment_variable() {
  local name=$1

  if [[ -z ${!name-} ]]; then
    die "Required environment variable ${name} is not configured."
  fi
  if [[ ${!name} == *$'\n'* || ${!name} == *$'\r'* ]]; then
    die "Required environment variable ${name} contains a line break."
  fi
}

valid_seqera_id() {
  local value=$1

  [[ $value =~ ^[A-Za-z0-9_-]+$ && ${#value} -le 128 ]]
}

api_json_request() {
  local label=$1
  local method=$2
  local path=$3
  local request_body=$4
  local response_file=$5
  local request_timeout=${6:-120}
  local -a curl_arguments=(
    --disable
    --fail
    --silent
    --show-error
    --proto '=https'
    --tlsv1.2
    --connect-timeout 30
    --max-time "$request_timeout"
    --request "$method"
    --url "${seqera_api_base}${path}"
    --header "@${auth_header_file}"
    --header 'Accept: application/json'
    --header 'Accept-Version: 1'
    --output "$response_file"
  )

  if [[ -n $request_body ]]; then
    curl_arguments+=(
      --header 'Content-Type: application/json'
      --data-binary "@${request_body}"
    )
  fi

  if ! curl "${curl_arguments[@]}"; then
    error "Seqera API ${label} request failed."
    return 1
  fi
  if ! jq -e . "$response_file" >/dev/null 2>&1; then
    error "Seqera API ${label} returned malformed JSON."
    return 1
  fi
}

cleanup_api_request() {
  local label=$1
  local method=$2
  local path=$3
  local request_body=${4-}
  local -a curl_arguments=(
    --disable
    --fail
    --silent
    --show-error
    --proto '=https'
    --tlsv1.2
    --connect-timeout 10
    --max-time 30
    --request "$method"
    --url "${seqera_api_base}${path}"
    --header "@${auth_header_file}"
    --header 'Accept: application/json'
    --header 'Accept-Version: 1'
    --output /dev/null
  )

  if [[ -n $request_body ]]; then
    curl_arguments+=(
      --header 'Content-Type: application/json'
      --data "$request_body"
    )
  fi

  if ! curl "${curl_arguments[@]}"; then
    error "Seqera API ${label} cleanup request failed."
    return 1
  fi
}

cleanup() {
  local original_status=$?
  local cleanup_status=0
  local recovered_workflow_id=''
  local recovery_response=''

  trap - EXIT
  trap '' HUP INT TERM
  set +e

  if [[ -z $workflow_id && -n $action_id && $launch_attempted -eq 1 ]]; then
    recovery_response="${scratch_dir}/describe-action.response.json"
    if api_json_request \
      'describe action during cleanup' \
      GET \
      "/actions/${action_id}?workspaceId=${TOWER_WORKSPACE_ID}" \
      '' \
      "$recovery_response" \
      20; then
      recovered_workflow_id=$(
        jq -er '.action.event.workflowId | select(type == "string")' \
          "$recovery_response" 2>/dev/null
      )
      if valid_seqera_id "$recovered_workflow_id"; then
        workflow_id=$recovered_workflow_id
        printf 'Recovered the launched workflow identifier during cleanup.\n'
      fi
    fi
  fi

  if [[ -n $workflow_id && $workflow_terminal -eq 0 ]]; then
    if cleanup_api_request \
      'cancel workflow' \
      POST \
      "/workflow/${workflow_id}/cancel?workspaceId=${TOWER_WORKSPACE_ID}" \
      '{}'; then
      printf 'Requested cancellation of incomplete Seqera workflow.\n'
    else
      cleanup_status=1
    fi
  fi

  if [[ -n $action_id ]]; then
    if cleanup_api_request \
      'delete action' \
      DELETE \
      "/actions/${action_id}?workspaceId=${TOWER_WORKSPACE_ID}"; then
      printf 'Deleted temporary Seqera action.\n'
    else
      cleanup_status=1
    fi
  fi

  if [[ -n $scratch_dir ]]; then
    rm -rf -- "$scratch_dir"
  fi

  if [[ $original_status -eq 0 && $cleanup_status -ne 0 ]]; then
    exit "$cleanup_status"
  fi
  exit "$original_status"
}

interrupted() {
  local signal_name=$1

  error "Seqera release gate interrupted by ${signal_name}."
  exit 1
}

trap cleanup EXIT
trap 'interrupted HUP' HUP
trap 'interrupted INT' INT
trap 'interrupted TERM' TERM

for required_name in \
  TOWER_AUTH_TOKEN \
  TOWER_WORK_DIR \
  PIPELINE_OUTDIR_ROOT \
  TOWER_WORKSPACE_ID \
  TOWER_COMPUTE_ENV_ID \
  CIRCLE_SHA1 \
  CIRCLE_BRANCH \
  CIRCLE_BUILD_NUM; do
  require_environment_variable "$required_name"
done

[[ $TOWER_WORKSPACE_ID =~ ^[0-9]+$ ]] || \
  die 'TOWER_WORKSPACE_ID must be a numeric workspace identifier.'
valid_seqera_id "$TOWER_COMPUTE_ENV_ID" || \
  die 'TOWER_COMPUTE_ENV_ID is not a valid Seqera identifier.'
[[ $CIRCLE_SHA1 =~ ^[0-9a-fA-F]{40}$ ]] || \
  die 'CIRCLE_SHA1 must be a full 40-character Git commit hash.'
[[ $CIRCLE_BUILD_NUM =~ ^[0-9]+$ ]] || \
  die 'CIRCLE_BUILD_NUM must be numeric.'
[[ ${#CIRCLE_BRANCH} -le 255 ]] || \
  die 'CIRCLE_BRANCH exceeds the supported length.'

readonly poll_interval_seconds=${SEQERA_POLL_INTERVAL_SECONDS:-60}
# Heavy validated fixtures have taken over three hours. Allow 4.5 hours while
# retaining a 30-minute cleanup margin below CircleCI Scale's five-hour job cap.
readonly poll_timeout_seconds=${SEQERA_POLL_TIMEOUT_SECONDS:-16200}
[[ $poll_interval_seconds =~ ^[0-9]+$ ]] || \
  die 'SEQERA_POLL_INTERVAL_SECONDS must be a non-negative integer.'
[[ $poll_timeout_seconds =~ ^[0-9]+$ ]] || \
  die 'SEQERA_POLL_TIMEOUT_SECONDS must be a non-negative integer.'

checked_out_sha=$(git rev-parse --verify HEAD 2>/dev/null) || \
  die 'Unable to resolve the checked-out Git commit.'
if [[ ${checked_out_sha,,} != "${CIRCLE_SHA1,,}" ]]; then
  die 'CIRCLE_SHA1 does not match the checked-out Git commit.'
fi

umask 077
scratch_dir=$(mktemp -d)
auth_header_file="${scratch_dir}/authorization.header"
printf 'Authorization: Bearer %s\n' "$TOWER_AUTH_TOKEN" > "$auth_header_file"

readonly action_response="${scratch_dir}/create-action.response.json"
readonly launch_response="${scratch_dir}/launch-action.response.json"
readonly workflow_response="${scratch_dir}/workflow.response.json"
readonly action_payload="${scratch_dir}/create-action.request.json"
readonly launch_payload="${scratch_dir}/launch-action.request.json"
readonly run_name="CI_${CIRCLE_SHA1:0:7}_${CIRCLE_BUILD_NUM}"
readonly output_directory="${PIPELINE_OUTDIR_ROOT%/}/circleci/${run_name}"

if ! jq -n \
  --arg name "$run_name" \
  --arg compute_environment_id "$TOWER_COMPUTE_ENV_ID" \
  --arg pipeline "$pipeline_url" \
  --arg work_directory "$TOWER_WORK_DIR" \
  --arg revision "$CIRCLE_BRANCH" \
  --arg commit_id "$CIRCLE_SHA1" \
  --arg nextflow_version "$nextflow_version" \
  '{
    name: $name,
    source: "tower",
    launch: {
      configProfiles: ["test"],
      computeEnvId: $compute_environment_id,
      pipeline: $pipeline,
      workDir: $work_directory,
      revision: $revision,
      commitId: $commit_id,
      pullLatest: false,
      preRunScript: ("export NXF_VER=" + $nextflow_version)
    }
  }' > "$action_payload"; then
  die 'Unable to construct the Seqera action request.'
fi

api_json_request \
  'create action' \
  POST \
  "/actions?workspaceId=${TOWER_WORKSPACE_ID}" \
  "$action_payload" \
  "$action_response" || exit 1

if ! candidate_action_id=$(
  jq -er '.actionId | select(type == "string")' "$action_response" 2>/dev/null
); then
  die 'Seqera create-action response did not contain a string actionId.'
fi
valid_seqera_id "$candidate_action_id" || \
  die 'Seqera create-action response contained an invalid actionId.'
action_id=$candidate_action_id

if ! jq -n \
  --arg output_directory "$output_directory" \
  '{params: {outdir: $output_directory}}' > "$launch_payload"; then
  die 'Unable to construct the Seqera launch request.'
fi

launch_attempted=1
api_json_request \
  'launch action' \
  POST \
  "/actions/${action_id}/launch?workspaceId=${TOWER_WORKSPACE_ID}" \
  "$launch_payload" \
  "$launch_response" || exit 1

if ! candidate_workflow_id=$(
  jq -er '.workflowId | select(type == "string")' "$launch_response" 2>/dev/null
); then
  die 'Seqera launch-action response did not contain a string workflowId.'
fi
valid_seqera_id "$candidate_workflow_id" || \
  die 'Seqera launch-action response contained an invalid workflowId.'
workflow_id=$candidate_workflow_id

printf 'Seqera workflow submitted for the exact CircleCI commit.\n'
SECONDS=0

while true; do
  api_json_request \
    'describe workflow' \
    GET \
    "/workflow/${workflow_id}?workspaceId=${TOWER_WORKSPACE_ID}" \
    '' \
    "$workflow_response" || exit 1

  if ! workflow_status=$(
    jq -er '.workflow.status | select(type == "string")' \
      "$workflow_response" 2>/dev/null
  ); then
    die 'Seqera describe-workflow response did not contain a string workflow status.'
  fi
  if ! observed_commit_id=$(
    jq -er '(.workflow.commitId // "") | select(type == "string")' \
      "$workflow_response" 2>/dev/null
  ); then
    die 'Seqera describe-workflow response contained a malformed commitId.'
  fi
  case "$workflow_status" in
    SUCCEEDED|FAILED|CANCELLED) workflow_terminal=1 ;;
  esac
  if [[ -n $observed_commit_id && ${observed_commit_id,,} != "${CIRCLE_SHA1,,}" ]]; then
    die 'Seqera reported a workflow commit that differs from CIRCLE_SHA1.'
  fi

  printf 'Seqera workflow status: %s (elapsed %ss).\n' "$workflow_status" "$SECONDS"
  case "$workflow_status" in
    SUCCEEDED)
      if [[ -z $observed_commit_id ]]; then
        die 'Succeeded Seqera workflow did not report its executed commitId.'
      fi
      printf 'Seqera workflow succeeded at the exact CircleCI commit.\n'
      break
      ;;
    FAILED)
      die 'Seqera workflow failed.'
      ;;
    CANCELLED)
      die 'Seqera workflow was cancelled.'
      ;;
    UNKNOWN)
      die 'Seqera workflow entered UNKNOWN status.'
      ;;
    SUBMITTED|RUNNING)
      if (( SECONDS >= poll_timeout_seconds )); then
        die "Seqera workflow exceeded the ${poll_timeout_seconds}s polling deadline."
      fi
      ;;
    *)
      die 'Seqera workflow returned an unexpected status.'
      ;;
  esac

  if (( poll_interval_seconds > 0 )); then
    sleep "$poll_interval_seconds"
  fi
done
