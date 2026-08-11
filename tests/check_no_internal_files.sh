#!/usr/bin/env bash
#
# Guard for a PUBLIC repository.
#
# Anything committed here is world-readable on push and remains readable in the git
# history after deletion. This script fails the build on the mechanically-detectable
# ways internal material has previously leaked in:
#
#   1. agent/AI working files and ad-hoc internal test data committed by path
#   2. credentials committed by content
#   3. beast-local filesystem paths committed by content
#
# Judgement calls that cannot be mechanised (internal S3 bucket names, customer/site
# run identifiers, colleague names) are covered by the checklist in CLAUDE.md, not here:
# the pipeline legitimately references csg-* buckets in .circleci/config.yml, so a
# content check on those would be noise rather than signal.
#
# Run from the repo root:  tests/check_no_internal_files.sh

set -euo pipefail
cd "$(git rev-parse --show-toplevel)"

SELF="tests/check_no_internal_files.sh"
status=0

fail() { printf '\nFAIL: %s\n' "$1"; status=1; }

# --- 1. Forbidden paths -------------------------------------------------------------
# Agent working files (planning docs, handoffs, campaign logs) and ad-hoc sample sheets
# pointing at internal data. These are gitignored, but .gitignore does not apply to files
# that are already tracked, so check the index directly.
forbidden_paths=$(git ls-files -- \
  'plans/**' 'plans' \
  'docker/agent/**' \
  '*_test_data/**' \
  'HANDOFF.md' '**/HANDOFF.md' \
  'mixed_testing.csv' \
  'input_csv/badsamples.csv' || true)

if [ -n "$forbidden_paths" ]; then
  fail "internal-only paths are tracked in this public repo:"
  printf '  %s\n' $forbidden_paths
  printf '\n  Agent working files belong in the session scratchpad or a private location.\n'
  printf '  Internal sample sheets belong in\n'
  printf '  s3://csg-reference/internal_nf_tests_data/external_validation/.\n'
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

if hits=$(git grep -n -I -E "$cred_re" -- . ":!$SELF" 2>/dev/null); then
  fail "possible credential committed:"
  printf '%s\n' "$hits" | sed 's/^/  /'
  printf '\n  If this is real, ROTATE IT FIRST -- removing the commit does not unpublish it.\n'
fi

# --- 3. Local machine paths by content ----------------------------------------------
# beast-local paths are meaningless to users and expose internal layout. CLAUDE.md and
# .gitignore document the rule and necessarily quote the patterns, so exclude them.
if hits=$(git grep -n -I -E '/nssd2/|agent-secrets|\.auth_tokens' \
            -- . ":!$SELF" ':!CLAUDE.md' ':!.gitignore' 2>/dev/null); then
  fail "beast-local path or secrets location committed:"
  printf '%s\n' "$hits" | sed 's/^/  /'
fi

if [ "$status" -eq 0 ]; then
  echo "OK: no internal-only paths, credentials or local paths tracked."
fi
exit "$status"
