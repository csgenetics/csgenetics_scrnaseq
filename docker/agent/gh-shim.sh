#!/usr/bin/env bash
# Thin wrapper around /usr/bin/gh (the real GitHub CLI).
#
# `gh` reads GH_TOKEN from the environment once at startup and does not
# re-read it.  Since we want a fresh installation access token on every
# invocation, this shim mints one via gh-credential-helper.sh and exports
# it before exec'ing the real gh.
#
# This file is installed at /usr/local/bin/gh by the Dockerfile, which
# shadows /usr/bin/gh in PATH so any `gh ...` call in the container hits
# the shim first.
set -euo pipefail

token="$(/usr/local/bin/gh-credential-helper.sh)"
exec env "GH_TOKEN=$token" /usr/bin/gh "$@"
