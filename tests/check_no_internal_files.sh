#!/usr/bin/env bash
#
# Guard for a PUBLIC repository.
#
# Anything committed here is world-readable on push. This script checks the Git index,
# which is the proposed commit for a local run and the checked-out commit tip in CI. It
# detects mechanically recognisable internal material staged by path or content:
#
#   1. agent/AI working files and ad-hoc internal test data
#   2. credentials
#   3. beast-local filesystem paths
#   4. private S3 locations expressed as s3/s3a/s3n URIs or standard AWS HTTPS endpoints
#
# Customer/site run identifiers and colleague names still require human review. Literal
# S3 buckets must identify themselves as public (with a delimited "public" component) or
# use the example convention. Private runtime locations must come from the environment.
# This guard does not inspect all history; already-published history cannot be retracted,
# so any historical credential exposure still requires rotation.
#
# Run from the repo root:  tests/check_no_internal_files.sh

set -euo pipefail
cd "$(git rev-parse --show-toplevel)"

status=0

fail() { printf '\nFAIL: %s\n' "$1"; status=1; }

print_match_locations() {
  local match remainder

  while IFS= read -r match; do
    [ -z "$match" ] && continue
    remainder=${match#*:}
    printf '  %s:%s\n' "${match%%:*}" "${remainder%%:*}"
  done
}

# git grep returns 1 for the expected "no matches" case and >1 for a real error. Keep
# those outcomes distinct so a broken scan cannot silently make the guard pass.
git_grep_index() {
  local output grep_status

  set +e
  output=$(git grep --cached "$@" 2>&1)
  grep_status=$?
  set -e
  if [ "$grep_status" -gt 1 ]; then
    printf 'ERROR: unable to scan the Git index:\n%s\n' "$output" >&2
    return "$grep_status"
  fi
  printf '%s' "$output"
}

# --- 1. Forbidden paths -------------------------------------------------------------
# Ignore rules do not apply to already-staged files, so check the index directly.
forbidden_paths=$(git ls-files -- \
  'plans/**' 'plans' \
  'docker/agent/**' \
  '*_test_data/**' \
  'HANDOFF.md' '**/HANDOFF.md' \
  'mixed_testing.csv' \
  'input_csv/badsamples.csv')

if [ -n "$forbidden_paths" ]; then
  fail "internal-only paths are staged in this public repo:"
  printf '  %s\n' $forbidden_paths
  printf '\n  Agent working files belong in the session scratchpad or a private location.\n'
  printf '  Internal sample sheets belong in access-controlled storage outside this repo.\n'
fi

# --- 2. Credentials by content ------------------------------------------------------
# Deliberately narrow: only patterns that are credentials and nothing else, so a hit is
# never a false positive.
cred_re='(AKIA[0-9A-Z]{16})'
cred_re+='|(ghp_[A-Za-z0-9]{20,})'
cred_re+='|(github_pat_[A-Za-z0-9_]{20,})'
cred_re+='|(sk-ant-[A-Za-z0-9-]{20,})'
cred_re+='|(-----BEGIN [A-Z ]*PRIVATE KEY)'
cred_re+='|(eyJ[A-Za-z0-9_-]{20,}\.[A-Za-z0-9_-]{20,}\.)'

hits=$(git_grep_index -n -I -E "$cred_re" -- .)
if [ -n "$hits" ]; then
  fail "possible credential staged:"
  print_match_locations <<< "$hits"
  printf '\n  If this is real, ROTATE IT FIRST -- removing the commit does not unpublish it.\n'
fi

# --- 3. Private S3 locations by content ---------------------------------------------
# Scope: AWS CLI/Hadoop URI schemes plus standard virtual-hosted and path-style HTTPS
# endpoints under amazonaws.com (including regional and China suffixes). Other cloud
# providers and custom S3-compatible domains remain human-review concerns.
is_public_or_example_bucket() {
  local bucket=$1
  local padded=".${bucket}."

  [[ "$padded" == *[._-]public[._-]* ]] ||
    [[ "$bucket" == example ]] ||
    [[ "$bucket" == example-* ]] ||
    [[ "$bucket" == *.example ]] ||
    [[ "$bucket" == *.example.* ]]
}

record_private_bucket() {
  local match=$1
  local bucket=${2,,}
  local remainder

  if ! is_public_or_example_bucket "$bucket"; then
    remainder=${match#*:}
    private_s3_hits+="${match%%:*}:${remainder%%:*}"$'\n'
  fi
}

private_s3_hits=''

s3_uri_re='s3[an]?://[A-Za-z0-9][A-Za-z0-9.-]{1,61}[A-Za-z0-9]'
uri_matches=$(git_grep_index -n -I -i -o -E "$s3_uri_re" -- .)
while IFS= read -r match; do
  [ -z "$match" ] && continue
  endpoint=${match#*:*:}
  record_private_bucket "$match" "${endpoint#*://}"
done < <(printf '%s\n' "$uri_matches")

virtual_host_re='https://[A-Za-z0-9][A-Za-z0-9.-]{1,61}[A-Za-z0-9]\.s3([.-][A-Za-z0-9-]+)*\.amazonaws\.com(\.cn)?'
virtual_host_matches=$(git_grep_index -n -I -i -o -E "$virtual_host_re" -- .)
while IFS= read -r match; do
  [ -z "$match" ] && continue
  endpoint=${match#*:*:}
  host=${endpoint#*://}
  host=${host,,}
  host_core=${host%.amazonaws.com.cn}
  host_core=${host_core%.amazonaws.com}
  service_suffix=${host_core##*.s3}
  bucket=${host_core%".s3${service_suffix}"}
  record_private_bucket "$match" "$bucket"
done < <(printf '%s\n' "$virtual_host_matches")

path_style_re='https://s3([.-][A-Za-z0-9-]+)*\.amazonaws\.com(\.cn)?/[A-Za-z0-9][A-Za-z0-9.-]{1,61}[A-Za-z0-9]'
path_style_matches=$(git_grep_index -n -I -i -o -E "$path_style_re" -- .)
while IFS= read -r match; do
  [ -z "$match" ] && continue
  endpoint=${match#*:*:}
  record_private_bucket "$match" "${endpoint##*/}"
done < <(printf '%s\n' "$path_style_matches")

if [ -n "$private_s3_hits" ]; then
  fail "private S3 endpoint literal staged:"
  printf '%s' "$private_s3_hits" | sed 's/^/  /'
  printf '\n  Inject private storage locations through the CI/runtime environment.\n'
fi

# --- 4. Local machine paths by content ----------------------------------------------
# CLAUDE.md and .gitignore document the rule and necessarily quote these patterns, so
# exclude those two policy files. Build the expression in pieces so this script itself
# remains scan-safe.
local_path_re='/nssd'
local_path_re+='2/|agent-'
local_path_re+='secrets|\.auth_'
local_path_re+='tokens'
hits=$(git_grep_index -n -I -E "$local_path_re" \
         -- . ':!CLAUDE.md' ':!.gitignore')
if [ -n "$hits" ]; then
  fail "beast-local path or secrets location staged:"
  print_match_locations <<< "$hits"
fi

if [ "$status" -eq 0 ]; then
  echo "OK: no internal-only paths, credentials, private S3 endpoints or local paths staged."
fi
exit "$status"
