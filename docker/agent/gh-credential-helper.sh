#!/usr/bin/env bash
# Mints a fresh GitHub App installation access token using the read-only
# private key mounted at /run/secrets/github-app.pem, and prints it in
# git credential helper format (on "get") or as a bare token (when called
# by the gh shim with no arguments).
#
# "store" and "erase" invocations (also part of the git credential helper
# protocol) are no-ops: this helper mints fresh tokens per call, so there
# is nothing to persist or clear. Short-circuiting avoids burning an
# install-token API call on every successful git push.
#
# Required env vars:
#   GITHUB_APP_ID
#   GITHUB_APP_INSTALLATION_ID
#
# Used by:
#   - git: as `credential.https://github.com.helper`
#   - gh: via the /usr/local/bin/gh shim
set -euo pipefail

# Git credential protocol: short-circuit on store/erase before doing any
# work. This helper is stateless — nothing to store or erase — so each
# such call would otherwise burn an install-token API call for no reason.
case "${1:-}" in
  store|erase) exit 0 ;;
esac

: "${GITHUB_APP_ID:?GITHUB_APP_ID is required}"
: "${GITHUB_APP_INSTALLATION_ID:?GITHUB_APP_INSTALLATION_ID is required}"
PEM="/run/secrets/github-app.pem"

b64url() { openssl base64 -A | tr -d '=' | tr '/+' '_-'; }

now=$(date +%s); iat=$((now - 60)); exp=$((now + 540))
header='{"typ":"JWT","alg":"RS256"}'
payload=$(printf '{"iat":%d,"exp":%d,"iss":"%s"}' "$iat" "$exp" "$GITHUB_APP_ID")
unsigned="$(printf '%s' "$header" | b64url).$(printf '%s' "$payload" | b64url)"
sig=$(printf '%s' "$unsigned" | openssl dgst -sha256 -sign "$PEM" -binary | b64url)
jwt="${unsigned}.${sig}"

resp=$(curl -sSL -X POST \
  -H "Authorization: Bearer $jwt" \
  -H "Accept: application/vnd.github+json" \
  -H "X-GitHub-Api-Version: 2022-11-28" \
  "https://api.github.com/app/installations/${GITHUB_APP_INSTALLATION_ID}/access_tokens")

token="$(echo "$resp" | jq -r '.token // empty')"
if [ -z "$token" ]; then
  echo "ERROR: failed to mint installation token" >&2
  echo "$resp" >&2
  exit 1
fi

if [ "${1:-}" = "get" ]; then
  printf 'username=x-access-token\npassword=%s\n' "$token"
else
  printf '%s\n' "$token"
fi
