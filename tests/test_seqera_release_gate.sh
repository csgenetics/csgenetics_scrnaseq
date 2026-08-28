#!/usr/bin/env bash

# Mocked contract tests for the CircleCI -> Seqera release gate. The test file
# doubles as a fake curl executable when invoked through the temporary `curl`
# symlink created below, so no network request can escape these tests.

set -euo pipefail

mock_curl() {
  local call_number=0
  local method=''
  local url=''
  local output_file=''
  local data_argument=''
  local header_argument=''
  local exit_status=0
  local saw_fail=0
  local saw_disable=0
  local saw_silent=0
  local saw_show_error=0
  local saw_https_only=0

  if [[ -f $MOCK_CURL_STATE/call-count ]]; then
    read -r call_number < "$MOCK_CURL_STATE/call-count"
  fi
  call_number=$((call_number + 1))
  printf '%s\n' "$call_number" > "$MOCK_CURL_STATE/call-count"

  while (( $# > 0 )); do
    case "$1" in
      --request)
        method=$2
        shift 2
        ;;
      --url)
        url=$2
        shift 2
        ;;
      --output)
        output_file=$2
        shift 2
        ;;
      --data|--data-binary)
        data_argument=$2
        shift 2
        ;;
      --header)
        header_argument=$2
        if [[ $header_argument == Authorization:* ]]; then
          : > "$MOCK_CURL_STATE/literal-authorization-header"
        elif [[ $header_argument == @* ]]; then
          header_file=${header_argument#@}
          if [[ ! -f $header_file ]] || \
            [[ $(<"$header_file") != "Authorization: Bearer ${MOCK_CURL_EXPECTED_TOKEN}" ]]; then
            : > "$MOCK_CURL_STATE/invalid-authorization-header"
          fi
        fi
        shift 2
        ;;
      --fail)
        saw_fail=1
        shift
        ;;
      --disable)
        saw_disable=1
        shift
        ;;
      --silent)
        saw_silent=1
        shift
        ;;
      --show-error)
        saw_show_error=1
        shift
        ;;
      --proto)
        [[ $2 == '=https' ]] && saw_https_only=1
        shift 2
        ;;
      *)
        shift
        ;;
    esac
  done

  printf '%s\t%s\t%s\n' "$call_number" "$method" "$url" \
    >> "$MOCK_CURL_STATE/calls.tsv"
  if [[ $saw_disable -ne 1 || $saw_fail -ne 1 || $saw_silent -ne 1 || \
        $saw_show_error -ne 1 || \
        $saw_https_only -ne 1 ]]; then
    : > "$MOCK_CURL_STATE/unsafe-curl-options"
  fi

  if [[ -n $data_argument ]]; then
    if [[ $data_argument == @* ]]; then
      cp -- "${data_argument#@}" "$MOCK_CURL_STATE/request-body.${call_number}.json"
    else
      printf '%s' "$data_argument" > "$MOCK_CURL_STATE/request-body.${call_number}.json"
    fi
  fi

  if [[ -f $MOCK_CURL_STATE/response.${call_number} ]]; then
    cp -- "$MOCK_CURL_STATE/response.${call_number}" "$output_file"
  else
    : > "$output_file"
  fi
  if [[ -f $MOCK_CURL_STATE/http-error.${call_number} && $saw_fail -eq 1 ]]; then
    exit_status=22
  elif [[ -f $MOCK_CURL_STATE/exit.${call_number} ]]; then
    read -r exit_status < "$MOCK_CURL_STATE/exit.${call_number}"
  fi
  return "$exit_status"
}

if [[ ${0##*/} == curl ]]; then
  set +e
  mock_curl "$@"
  exit $?
fi

repo_root=$(git rev-parse --show-toplevel)
scratch_root=$(mktemp -d)
trap 'rm -rf -- "$scratch_root"' EXIT
mock_bin="${scratch_root}/mock-bin"
mkdir -p "$mock_bin"
ln -s "$repo_root/tests/test_seqera_release_gate.sh" "$mock_bin/curl"
ln -s /bin/true "$mock_bin/sleep"

readonly expected_token='mock-seqera-token-sentinel'
readonly private_work_directory='mock://private-work-value-"quoted"'
readonly private_output_root='mock://private-output-value-"quoted"'
readonly workspace_id='12345678901234'
readonly compute_environment_id='compute_ABC-123'
expected_sha=$(git -C "$repo_root" rev-parse HEAD)
readonly expected_sha
readonly expected_branch='release/mock-branch'

grep -Fq 'SEQERA_POLL_TIMEOUT_SECONDS:-16200' \
  "$repo_root/.circleci/run_seqera_release_gate.sh" || {
    printf 'FAIL: production Seqera polling deadline is not 4.5 hours\n' >&2
    exit 1
  }

case_state=''
case_output=''
case_status=0
case_timeout=60
case_poll_interval=1
case_sha=$expected_sha
case_circle_sha=$expected_sha
case_circle_branch=$expected_branch
case_revision=$expected_branch
include_token=1

new_case() {
  local name=$1

  case_state="${scratch_root}/${name}"
  mkdir -p "$case_state"
  case_timeout=60
  case_poll_interval=1
  case_sha=$expected_sha
  case_circle_sha=$expected_sha
  case_circle_branch=$expected_branch
  case_revision=$expected_branch
  include_token=1
}

response() {
  local call_number=$1
  local body=$2

  printf '%s' "$body" > "$case_state/response.${call_number}"
}

http_failure() {
  local call_number=$1

  : > "$case_state/http-error.${call_number}"
}

run_case() {
  local -a environment=(
    "PATH=${mock_bin}:${PATH}"
    "MOCK_CURL_STATE=${case_state}"
    "MOCK_CURL_EXPECTED_TOKEN=${expected_token}"
    "TOWER_WORK_DIR=${private_work_directory}"
    "PIPELINE_OUTDIR_ROOT=${private_output_root}"
    "TOWER_WORKSPACE_ID=${workspace_id}"
    "TOWER_COMPUTE_ENV_ID=${compute_environment_id}"
    "CIRCLE_SHA1=${case_circle_sha}"
    "CIRCLE_BRANCH=${case_circle_branch}"
    "SEQERA_COMMIT_ID=${case_sha}"
    "SEQERA_REVISION=${case_revision}"
    'CIRCLE_BUILD_NUM=4242'
    "SEQERA_POLL_INTERVAL_SECONDS=${case_poll_interval}"
    "SEQERA_POLL_TIMEOUT_SECONDS=${case_timeout}"
  )

  if [[ $include_token -eq 1 ]]; then
    environment+=("TOWER_AUTH_TOKEN=${expected_token}")
  fi

  set +e
  case_output=$(
    cd "$repo_root" && \
      env -i "${environment[@]}" bash .circleci/run_seqera_release_gate.sh 2>&1
  )
  case_status=$?
  set -e
}

fail_test() {
  local message=$1

  printf 'FAIL: %s\n' "$message" >&2
  if [[ -n $case_output ]]; then
    printf '%s\n' "$case_output" >&2
  fi
  exit 1
}

expect_success() {
  local name=$1

  if [[ $case_status -ne 0 ]]; then
    fail_test "${name} should have succeeded"
  fi
}

expect_failure() {
  local name=$1
  local diagnostic=$2

  if [[ $case_status -eq 0 ]]; then
    fail_test "${name} should have failed"
  fi
  if [[ $case_output != *"$diagnostic"* ]]; then
    fail_test "${name} lacked its expected diagnostic: ${diagnostic}"
  fi
}

assert_private_values_hidden() {
  local name=$1

  for private_value in \
    "$expected_token" \
    "$private_work_directory" \
    "$private_output_root"; do
    if [[ $case_output == *"$private_value"* ]]; then
      fail_test "${name} printed a secret or private value"
    fi
  done
}

assert_header_was_file_backed() {
  local name=$1

  if [[ -e $case_state/literal-authorization-header || \
        -e $case_state/invalid-authorization-header ]]; then
    fail_test "${name} did not pass the authorization header through its private file"
  fi
}

assert_curl_was_fail_loud() {
  local name=$1

  if [[ -e $case_state/unsafe-curl-options ]]; then
    fail_test "${name} omitted a required fail-loud curl option"
  fi
}

assert_call() {
  local call_number=$1
  local method=$2
  local path=$3

  if ! awk -F '\t' \
    -v number="$call_number" \
    -v expected_method="$method" \
    -v expected_url="https://api.cloud.seqera.io${path}" \
    '$1 == number && $2 == expected_method && $3 == expected_url { found = 1 }
     END { exit(found ? 0 : 1) }' \
    "$case_state/calls.tsv"; then
    fail_test "missing mocked API call ${call_number}: ${method} ${path}"
  fi
}

# Success: payloads preserve special characters through jq, the launch is pinned
# to the exact commit, every active status is logged, and the action is deleted.
new_case success
response 1 '{"actionId":"action_ABC-123"}'
response 2 '{"workflowId":"workflow_XYZ-456"}'
response 3 "{\"workflow\":{\"status\":\"SUBMITTED\",\"commitId\":\"${expected_sha}\"}}"
response 4 "{\"workflow\":{\"status\":\"RUNNING\",\"commitId\":\"${expected_sha}\"}}"
response 5 "{\"workflow\":{\"status\":\"SUCCEEDED\",\"commitId\":\"${expected_sha}\"}}"
response 6 ''
run_case
expect_success success
assert_private_values_hidden success
assert_header_was_file_backed success
assert_curl_was_fail_loud success
assert_call 1 POST "/actions?workspaceId=${workspace_id}"
assert_call 2 POST "/actions/action_ABC-123/launch?workspaceId=${workspace_id}"
assert_call 3 GET "/workflow/workflow_XYZ-456?workspaceId=${workspace_id}"
assert_call 6 DELETE "/actions/action_ABC-123?workspaceId=${workspace_id}"
jq -e \
  --arg sha "$expected_sha" \
  --arg revision "$expected_branch" \
  --arg work_directory "$private_work_directory" \
  '.launch.commitId == $sha and
   .launch.revision == $revision and
   .launch.pullLatest == false and
   .launch.workDir == $work_directory and
   .launch.configProfiles == ["test"]' \
  "$case_state/request-body.1.json" >/dev/null || \
  fail_test 'success action payload did not preserve the exact launch contract'
jq -e \
  --arg expected "${private_output_root}/circleci/CI_${expected_sha:0:7}_4242" \
  '.params.outdir == $expected' \
  "$case_state/request-body.2.json" >/dev/null || \
  fail_test 'success launch payload did not preserve the private output root'
for expected_status in SUBMITTED RUNNING SUCCEEDED; do
  [[ $case_output == *"Seqera workflow status: ${expected_status}"* ]] || \
    fail_test "success output omitted the ${expected_status} heartbeat"
done

# A missing secret is rejected before any API request.
new_case missing_token
include_token=0
run_case
expect_failure missing_token 'Required environment variable TOWER_AUTH_TOKEN is not configured.'
assert_private_values_hidden missing_token
[[ ! -e $case_state/calls.tsv ]] || \
  fail_test 'missing-token case reached the API'

# A zero interval would hammer the API continuously for the whole deadline.
new_case zero_poll_interval
case_poll_interval=0
run_case
expect_failure zero_poll_interval \
  'SEQERA_POLL_INTERVAL_SECONDS must be a positive integer.'
assert_private_values_hidden zero_poll_interval
[[ ! -e $case_state/calls.tsv ]] || \
  fail_test 'zero-poll-interval case reached the API'

# A missing pin is rejected before any API request. CIRCLE_SHA1 is not the pin.
new_case missing_pin
case_sha=''
run_case
expect_failure missing_pin \
  'Required environment variable SEQERA_COMMIT_ID is not configured.'
assert_private_values_hidden missing_pin
[[ ! -e $case_state/calls.tsv ]] || \
  fail_test 'missing-pin case reached the API'

# HTTP failures fail loud without replaying the response body.
new_case create_http_failure
response 1 '{"message":"private-api-error-body-sentinel"}'
http_failure 1
run_case
expect_failure create_http_failure 'Seqera API create action request failed.'
assert_private_values_hidden create_http_failure
assert_curl_was_fail_loud create_http_failure
[[ $case_output != *'private-api-error-body-sentinel'* ]] || \
  fail_test 'HTTP failure replayed the API response body'

# A successful HTTP response with malformed JSON is never parsed as an ID.
new_case malformed_json
response 1 '{not-json'
run_case
expect_failure malformed_json 'Seqera API create action returned malformed JSON.'
assert_private_values_hidden malformed_json
[[ $(<"$case_state/call-count") == 1 ]] || \
  fail_test 'malformed JSON response triggered an API side effect'

# Malformed create responses cannot be mistaken for valid action identifiers.
new_case malformed_action_id
response 1 '{"actionId":42}'
run_case
expect_failure malformed_action_id \
  'Seqera create-action response did not contain a string actionId.'
assert_private_values_hidden malformed_action_id
[[ $(<"$case_state/call-count") == 1 ]] || \
  fail_test 'malformed action ID triggered an unsafe cleanup request'

# Once an action exists, a malformed launch response still deletes that action.
new_case malformed_workflow_id
response 1 '{"actionId":"action_ABC-123"}'
response 2 '{"workflowId":null}'
response 3 '{"action":{"event":null}}'
response 4 ''
run_case
expect_failure malformed_workflow_id \
  'Seqera launch-action response did not contain a string workflowId.'
assert_private_values_hidden malformed_workflow_id
assert_call 3 GET "/actions/action_ABC-123?workspaceId=${workspace_id}"
assert_call 4 DELETE "/actions/action_ABC-123?workspaceId=${workspace_id}"

# A launch HTTP error can mean the server accepted the request but its response
# was lost. Recover the workflow ID from the new action before cleanup.
new_case launch_http_failure
response 1 '{"actionId":"action_ABC-123"}'
response 2 '{"message":"launch-response-lost"}'
http_failure 2
response 3 '{"action":{"event":{"workflowId":"workflow_XYZ-456"}}}'
response 4 ''
response 5 ''
run_case
expect_failure launch_http_failure 'Seqera API launch action request failed.'
assert_private_values_hidden launch_http_failure
assert_curl_was_fail_loud launch_http_failure
assert_call 3 GET "/actions/action_ABC-123?workspaceId=${workspace_id}"
assert_call 4 POST "/workflow/workflow_XYZ-456/cancel?workspaceId=${workspace_id}"
assert_call 5 DELETE "/actions/action_ABC-123?workspaceId=${workspace_id}"

# A malformed status field is distinct from a valid but unsupported status and
# still receives complete cancellation/deletion cleanup.
new_case malformed_status
response 1 '{"actionId":"action_ABC-123"}'
response 2 '{"workflowId":"workflow_XYZ-456"}'
response 3 "{\"workflow\":{\"status\":null,\"commitId\":\"${expected_sha}\"}}"
response 4 ''
response 5 ''
run_case
expect_failure malformed_status \
  'Seqera describe-workflow response did not contain a string workflow status.'
assert_private_values_hidden malformed_status
assert_call 4 POST "/workflow/workflow_XYZ-456/cancel?workspaceId=${workspace_id}"
assert_call 5 DELETE "/actions/action_ABC-123?workspaceId=${workspace_id}"

# UNKNOWN and future/unexpected statuses are failures; an incomplete workflow is
# cancelled before its temporary action is deleted.
for status_case in UNKNOWN PAUSED; do
  new_case "status_${status_case}"
  response 1 '{"actionId":"action_ABC-123"}'
  response 2 '{"workflowId":"workflow_XYZ-456"}'
  response 3 "{\"workflow\":{\"status\":\"${status_case}\",\"commitId\":\"${expected_sha}\"}}"
  response 4 ''
  response 5 ''
  run_case
  if [[ $status_case == UNKNOWN ]]; then
    expect_failure "status_${status_case}" 'Seqera workflow entered UNKNOWN status.'
  else
    expect_failure "status_${status_case}" 'Seqera workflow returned an unexpected status.'
  fi
  assert_private_values_hidden "status_${status_case}"
  assert_call 4 POST "/workflow/workflow_XYZ-456/cancel?workspaceId=${workspace_id}"
  assert_call 5 DELETE "/actions/action_ABC-123?workspaceId=${workspace_id}"
  jq -e '. == {}' "$case_state/request-body.4.json" >/dev/null || \
    fail_test "status_${status_case} cancellation did not send the API's empty JSON body"
done

# A bounded polling deadline cancels and deletes rather than looping forever.
new_case timeout
case_timeout=0
response 1 '{"actionId":"action_ABC-123"}'
response 2 '{"workflowId":"workflow_XYZ-456"}'
response 3 "{\"workflow\":{\"status\":\"RUNNING\",\"commitId\":\"${expected_sha}\"}}"
response 4 ''
response 5 ''
run_case
expect_failure timeout 'Seqera workflow exceeded the 0s polling deadline.'
assert_private_values_hidden timeout
assert_call 4 POST "/workflow/workflow_XYZ-456/cancel?workspaceId=${workspace_id}"
assert_call 5 DELETE "/actions/action_ABC-123?workspaceId=${workspace_id}"

# A terminal workflow failure does not issue a meaningless cancellation, but its
# temporary action is still removed.
new_case workflow_failed
response 1 '{"actionId":"action_ABC-123"}'
response 2 '{"workflowId":"workflow_XYZ-456"}'
response 3 "{\"workflow\":{\"status\":\"FAILED\",\"commitId\":\"${expected_sha}\"}}"
response 4 ''
run_case
expect_failure workflow_failed 'Seqera workflow failed.'
assert_private_values_hidden workflow_failed
assert_call 4 DELETE "/actions/action_ABC-123?workspaceId=${workspace_id}"
[[ $(<"$case_state/call-count") == 4 ]] || \
  fail_test 'terminal failure issued an unnecessary cancellation request'

# Cleanup is part of the gate: deletion failure turns an otherwise successful run
# into a failed CircleCI job.
new_case delete_failure
response 1 '{"actionId":"action_ABC-123"}'
response 2 '{"workflowId":"workflow_XYZ-456"}'
response 3 "{\"workflow\":{\"status\":\"SUCCEEDED\",\"commitId\":\"${expected_sha}\"}}"
response 4 ''
http_failure 4
run_case
expect_failure delete_failure 'Seqera API delete action cleanup request failed.'
assert_private_values_hidden delete_failure

# Even a successful workflow must prove which commit actually ran.
new_case missing_executed_commit
response 1 '{"actionId":"action_ABC-123"}'
response 2 '{"workflowId":"workflow_XYZ-456"}'
response 3 '{"workflow":{"status":"SUCCEEDED","commitId":null}}'
response 4 ''
run_case
expect_failure missing_executed_commit \
  'Succeeded Seqera workflow did not report its executed commitId.'
assert_private_values_hidden missing_executed_commit
assert_call 4 DELETE "/actions/action_ABC-123?workspaceId=${workspace_id}"
[[ $(<"$case_state/call-count") == 4 ]] || \
  fail_test 'missing executed commit attempted to cancel a terminal workflow'

# A commit mismatch reported by Seqera fails the gate and cancels an active run.
new_case mismatched_executed_commit
response 1 '{"actionId":"action_ABC-123"}'
response 2 '{"workflowId":"workflow_XYZ-456"}'
response 3 '{"workflow":{"status":"RUNNING","commitId":"1111111111111111111111111111111111111111"}}'
response 4 ''
response 5 ''
run_case
expect_failure mismatched_executed_commit \
  'Seqera reported a workflow commit that differs from SEQERA_COMMIT_ID.'
assert_private_values_hidden mismatched_executed_commit
assert_call 4 POST "/workflow/workflow_XYZ-456/cancel?workspaceId=${workspace_id}"
assert_call 5 DELETE "/actions/action_ABC-123?workspaceId=${workspace_id}"

# On-demand pin: CIRCLE_SHA1 is devel HEAD and must not be the Seqera commit.
new_case on_demand_pin
readonly on_demand_sha='aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa'
case_sha=$on_demand_sha
case_circle_sha=$expected_sha
case_revision='feature/on-demand'
response 1 '{"actionId":"action_ABC-123"}'
response 2 '{"workflowId":"workflow_XYZ-456"}'
response 3 "{\"workflow\":{\"status\":\"SUCCEEDED\",\"commitId\":\"${on_demand_sha}\"}}"
response 4 ''
run_case
expect_success on_demand_pin
assert_private_values_hidden on_demand_pin
jq -e \
  --arg sha "$on_demand_sha" \
  --arg revision 'feature/on-demand' \
  '.launch.commitId == $sha and .launch.revision == $revision and .launch.pullLatest == false' \
  "$case_state/request-body.1.json" >/dev/null || \
  fail_test 'on-demand pin used CIRCLE_SHA1 or dropped pullLatest=false'

# Tag pipelines leave CIRCLE_BRANCH empty. The pin is SEQERA_REVISION.
new_case tag_revision
case_circle_branch=''
case_revision='2.0.0'
response 1 '{"actionId":"action_ABC-123"}'
response 2 '{"workflowId":"workflow_XYZ-456"}'
response 3 "{\"workflow\":{\"status\":\"SUCCEEDED\",\"commitId\":\"${expected_sha}\"}}"
response 4 ''
run_case
expect_success tag_revision
assert_private_values_hidden tag_revision
jq -e \
  --arg sha "$expected_sha" \
  --arg revision '2.0.0' \
  '.launch.commitId == $sha and .launch.revision == $revision' \
  "$case_state/request-body.1.json" >/dev/null || \
  fail_test 'tag revision was not used as the Seqera revision'

# --- config.yml contract: no checkout, one job, devel/version-tag filters, on-demand ---
config_yml="$repo_root/.circleci/config.yml"
launcher_file="$repo_root/.circleci/run_seqera_release_gate.sh"

[[ -f $config_yml ]] || fail_test 'config.yml is missing'
grep -Fq 'parameters:' "$config_yml" || \
  fail_test 'config.yml has no pipeline parameters'
grep -Fq 'launch_sha:' "$config_yml" || \
  fail_test 'config.yml is missing launch_sha'
grep -Fq 'seqera-on-demand:' "$config_yml" || \
  fail_test 'config.yml is missing the seqera-on-demand workflow'

[[ -f $repo_root/.github/workflows/branch.yml ]] || \
  fail_test 'branch.yml is missing'
grep -Fq 'github.head_ref == devel' "$repo_root/.github/workflows/branch.yml" || \
  grep -Fq 'GITHUB_HEAD_REF == devel' "$repo_root/.github/workflows/branch.yml" || \
  fail_test 'branch.yml does not require head ref devel'

python3 - "$config_yml" "$launcher_file" <<'PY' || fail_test 'config.yml credentialed-job contract failed'
import pathlib, re, sys

config = pathlib.Path(sys.argv[1]).read_text()
launcher = pathlib.Path(sys.argv[2]).read_text()

begin = "# BEGIN_SEQERA_RELEASE_GATE\n"
end = "# END_SEQERA_RELEASE_GATE"
if begin not in config or end not in config:
    raise SystemExit("missing BEGIN/END_SEQERA_RELEASE_GATE markers")
marked = config.split(begin, 1)[1].split(end, 1)[0]
# Strip the YAML block indent shared by every line of the marked region.
lines = marked.splitlines()
indents = [len(line) - len(line.lstrip(" ")) for line in lines if line.strip()]
if not indents:
    raise SystemExit("marked launcher body is empty")
strip = min(indents)
body = "\n".join(line[strip:] if len(line) >= strip else line for line in lines)
if body.startswith("\n"):
    body = body[1:]
if not body.endswith("\n"):
    body += "\n"
if body != launcher:
    raise SystemExit("marked launcher body is not byte-identical to run_seqera_release_gate.sh")

job = re.search(
    r"^  run-current-branch:\n(?:.*\n)*?^    steps:\n((?:^      .*\n)+)",
    config,
    re.M,
)
if not job:
    raise SystemExit("could not find run-current-branch steps")
steps = job.group(1)
if re.search(r"^\s*- checkout\s*$", steps, re.M):
    raise SystemExit("run-current-branch still checkouts")
if ".circleci/run_seqera_release_gate.sh" in steps:
    raise SystemExit("run-current-branch still executes the checkout script path")
if "SEQERA_COMMIT_ID" not in config.split(begin, 1)[0]:
    raise SystemExit("prelude does not set SEQERA_COMMIT_ID")

if not re.search(
    r'seqera-on-demand:[\s\S]*?equal:\s*\[\s*devel,\s*"<< pipeline.git.branch >>"\s*\]',
    config,
):
    raise SystemExit("seqera-on-demand when is not restricted to devel")
if re.search(
    r"seqera-on-demand:[\s\S]*?pipeline.git.tag",
    config,
):
    raise SystemExit("seqera-on-demand when still allows tags")

# Credentialed job entry in the ci workflow: devel branches + version tags.
ci = re.search(r"^  ci:\n([\s\S]*?)^  seqera-on-demand:", config, re.M)
if not ci:
    raise SystemExit("could not find ci workflow before seqera-on-demand")
ci_block = ci.group(1)
cred = re.search(
    r"- run-current-branch:\n([\s\S]*?)(?:^      - |\Z)",
    ci_block,
)
if not cred:
    raise SystemExit("ci workflow does not name run-current-branch")
cred_block = cred.group(1)
if "only:\n                - devel" not in cred_block and "only: devel" not in cred_block:
    if not re.search(r"branches:\n\s+only:\n\s+- devel", cred_block):
        raise SystemExit("ci run-current-branch is not filtered to devel")
if r"/^[0-9]+\.[0-9]+\.[0-9]+$/" not in cred_block:
    raise SystemExit("ci run-current-branch is not filtered to version tags")
if r"only: /.*/" in cred_block:
    raise SystemExit("ci run-current-branch still uses tags.only: /.*/")

# Every other job listed under ci must have a tags filter so version-tag
# pipelines can satisfy requires.
for match in re.finditer(r"^      - ([a-z0-9-]+):", ci_block, re.M):
    name = match.group(1)
    if name == "run-current-branch":
        continue
    # Slice from this job to the next job or end
    start = match.start()
    nxt = re.search(r"^      - [a-z0-9-]+:", ci_block[match.end():], re.M)
    chunk = ci_block[start: match.end() + (nxt.start() if nxt else len(ci_block))]
    if "tags:" not in chunk:
        raise SystemExit(f"ci job {name} has no tags filter")
PY

printf 'OK: CircleCI -> Seqera release-gate contract tests passed.\n'
