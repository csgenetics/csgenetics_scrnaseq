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
SECRETS_DIR="${SECRETS_DIR:?Set SECRETS_DIR to the directory containing agent secrets (aws/, claude-oauth-token, github/)}"
DOCKER_IMAGE="${DOCKER_IMAGE:-scrnaseq-agent:latest}"
# App ID and Installation ID are public identifiers of the GitHub App,
# not secrets. They are useful only in combination with the private key.
GITHUB_APP_ID="${GITHUB_APP_ID:-3409188}"
GITHUB_APP_INSTALLATION_ID="${GITHUB_APP_INSTALLATION_ID:-124688983}"

OAUTH_FILE="$SECRETS_DIR/claude-oauth-token"
APP_PEM="$SECRETS_DIR/github/scrnaseq-agent-app.pem"
AWS_DIR="$SECRETS_DIR/aws"

# ---- Validate every input we need ----
for f in "$OAUTH_FILE" "$APP_PEM"; do
  [ -f "$f" ] || { echo "ERROR: missing file $f" >&2; exit 1; }
done
[ -d "$AWS_DIR" ]   || { echo "ERROR: missing dir $AWS_DIR"   >&2; exit 1; }
[ -d "$REPO_PATH" ] || { echo "ERROR: missing dir $REPO_PATH" >&2; exit 1; }

# ---- Load Claude token from host file (never mounted into container) ----
CLAUDE_CODE_OAUTH_TOKEN="$(tr -d '\n\r' < "$OAUTH_FILE")"
export CLAUDE_CODE_OAUTH_TOKEN GITHUB_APP_ID GITHUB_APP_INSTALLATION_ID

# ---- Run container ----
exec docker run \
  --rm \
  --name "scrnaseq-agent-$$" \
  --user 1000:1000 \
  -v "$REPO_PATH:/scrnaseq_git_repo:rw" \
  -v "$AWS_DIR:/home/node/.aws:ro" \
  -v "$APP_PEM:/run/secrets/github-app.pem:ro" \
  -e CLAUDE_CODE_OAUTH_TOKEN \
  -e GITHUB_APP_ID \
  -e GITHUB_APP_INSTALLATION_ID \
  "$DOCKER_IMAGE" \
  "$@"
