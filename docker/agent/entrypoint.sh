#!/usr/bin/env bash
# Entrypoint for the csgenetics-scrnaseq-agent container.
#
# Required env vars (injected by docker run):
#   CLAUDE_CODE_OAUTH_TOKEN      Claude Max subscription OAuth token
#   GRATICULE_API_KEY            Per-repo bot key for `scrnaseq-agent` (GRT-385).
#                                Consumed by the MCP config's `${GRATICULE_API_KEY}`
#                                substitution in claude-state.json — every MCP
#                                request to Graticule flows through this key.
#   GITHUB_APP_ID                Numeric GitHub App ID
#   GITHUB_APP_INSTALLATION_ID   Installation ID for csgenetics/csgenetics_scrnaseq
#
# Required read-only mount:
#   /run/secrets/github-app.pem  RSA private key for the GitHub App
#
# Usage:
#   entrypoint.sh                       # runs a default diagnostic prompt
#   entrypoint.sh "Your prompt"         # runs with the given prompt
#   entrypoint.sh "Your prompt" --extra-claude-flag
set -euo pipefail

: "${CLAUDE_CODE_OAUTH_TOKEN:?CLAUDE_CODE_OAUTH_TOKEN is required}"
: "${GRATICULE_API_KEY:?GRATICULE_API_KEY is required}"
: "${GITHUB_APP_ID:?GITHUB_APP_ID is required}"
: "${GITHUB_APP_INSTALLATION_ID:?GITHUB_APP_INSTALLATION_ID is required}"
[ -r /run/secrets/github-app.pem ] || {
  echo "ERROR: /run/secrets/github-app.pem is not readable" >&2
  exit 1
}

# ---- Git identity (standard GitHub App bot format) ----
git config --global user.name  "csgenetics-scrnaseq-agent[bot]"
git config --global user.email "${GITHUB_APP_ID}+csgenetics-scrnaseq-agent[bot]@users.noreply.github.com"

# ---- Git credential helper, scoped to github.com only ----
git config --global credential.https://github.com.helper \
  "/usr/local/bin/gh-credential-helper.sh"

# ---- Diagnostic prompt default ----
if [ $# -eq 0 ]; then
  PROMPT='You are running inside the csgenetics-scrnaseq-agent container for diagnostics. Run each of these bash commands and report the outputs in a single JSON object with one key per command: `pwd`, `git branch --show-current`, `git log --oneline -1`, `git config user.name`, `git config user.email`, `whoami`, `id -u`, `ls /scrnaseq_git_repo | head -5`. Then exit. Do not modify any files.'
  EXTRA_ARGS=()
else
  PROMPT="$1"
  shift
  EXTRA_ARGS=("$@")
fi

cd /scrnaseq_git_repo

exec claude \
  -p "$PROMPT" \
  --model "claude-opus-4-6[1m]" \
  --output-format json \
  --dangerously-skip-permissions \
  "${EXTRA_ARGS[@]}"
