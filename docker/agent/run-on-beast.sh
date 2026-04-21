#!/usr/bin/env bash
# Host-side wrapper for running the scrnaseq-agent container.
#
# Validates that all required files exist, loads the Claude OAuth token,
# and launches `docker run` with the right mounts and env vars.
#
# Usage:
#   docker/agent/run-on-beast.sh                      # default diagnostic prompt
#   docker/agent/run-on-beast.sh "Your prompt here"
#   docker/agent/run-on-beast.sh "Your prompt" --extra-claude-flag
#
# REQUIRED env vars (no defaults — set these before running):
#   REPO_PATH                  absolute path to your csgenetics_scrnaseq checkout
#   SECRETS_DIR                directory containing aws/, claude-oauth-token,
#                              github/scrnaseq-agent-app.pem
#
# Optional env vars (sensible defaults provided):
#   DOCKER_IMAGE               image tag to run (default: scrnaseq-agent:latest)
#   GITHUB_APP_ID              numeric app ID (default: 3409188)
#   GITHUB_APP_INSTALLATION_ID installation ID (default: 124688983)
set -euo pipefail

REPO_PATH="${REPO_PATH:?Set REPO_PATH to the host path of your csgenetics_scrnaseq checkout}"
SECRETS_DIR="${SECRETS_DIR:?Set SECRETS_DIR to the directory containing agent secrets (aws/, claude-oauth-token, github/, graticule-api-key)}"
DOCKER_IMAGE="${DOCKER_IMAGE:-scrnaseq-agent:latest}"
# App ID and Installation ID are public identifiers of the GitHub App,
# not secrets. They are useful only in combination with the private key.
GITHUB_APP_ID="${GITHUB_APP_ID:-3409188}"
GITHUB_APP_INSTALLATION_ID="${GITHUB_APP_INSTALLATION_ID:-124688983}"
# Graticule base URL is a public identifier.
GRATICULE_BASE_URL="${GRATICULE_BASE_URL:-https://graticule.csgenetics.com}"

OAUTH_FILE="$SECRETS_DIR/claude-oauth-token"
# Per-repo Graticule bot API key (GRT-385). Authenticates as the
# `scrnaseq-agent` bot user, not the operator. The filename is a public
# identifier; only the file's contents are secret, and it never leaves
# SECRETS_DIR except as a bind-mount + env var into the container.
GRATICULE_API_KEY_FILE="$SECRETS_DIR/graticule-api-key"
APP_PEM="$SECRETS_DIR/github/scrnaseq-agent-app.pem"
AWS_DIR="$SECRETS_DIR/aws"

# ---- Validate every input we need ----
for f in "$OAUTH_FILE" "$GRATICULE_API_KEY_FILE" "$APP_PEM"; do
  [ -f "$f" ] || { echo "ERROR: missing file $f" >&2; exit 1; }
done
[ -d "$AWS_DIR" ]   || { echo "ERROR: missing dir $AWS_DIR"   >&2; exit 1; }
[ -d "$REPO_PATH" ] || { echo "ERROR: missing dir $REPO_PATH" >&2; exit 1; }

# ---- Load secrets from host files (never mounted into container) ----
CLAUDE_CODE_OAUTH_TOKEN="$(tr -d '\n\r' < "$OAUTH_FILE")"
GRATICULE_API_KEY="$(tr -d '\n\r' < "$GRATICULE_API_KEY_FILE")"
export CLAUDE_CODE_OAUTH_TOKEN GRATICULE_API_KEY GRATICULE_BASE_URL \
       GITHUB_APP_ID GITHUB_APP_INSTALLATION_ID

# ---- Run container ----
exec docker run \
  --rm \
  --name "scrnaseq-agent-$$" \
  --user 1000:1000 \
  -v "$REPO_PATH:/scrnaseq_git_repo:rw" \
  -v "$AWS_DIR:/home/node/.aws:ro" \
  -v "$APP_PEM:/run/secrets/github-app.pem:ro" \
  -e CLAUDE_CODE_OAUTH_TOKEN \
  -e GRATICULE_API_KEY \
  -e GRATICULE_BASE_URL \
  -e GITHUB_APP_ID \
  -e GITHUB_APP_INSTALLATION_ID \
  "$DOCKER_IMAGE" \
  "$@"
