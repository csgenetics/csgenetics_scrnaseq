#!/usr/bin/env bash

# Verify that the exact STAR container configured for production can execute
# the repository parser added to the STAR process. This complements the Conda
# dependency contract with a real container-path invocation.

set -euo pipefail

repo_root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
image_config="${repo_root}/conf/images.config"
fixture="${repo_root}/tests/fixtures/star/Log.final.out"

star_image=$(
  sed -n "s/.*withName: star[[:space:]]*{[[:space:]]*container = '\([^']*\)'.*/\1/p" \
    "$image_config"
)
[[ -n $star_image ]] || {
  printf 'ERROR: unable to resolve the configured STAR container image\n' >&2
  exit 1
}
[[ -f $fixture ]] || {
  printf 'ERROR: STAR runtime fixture is missing\n' >&2
  exit 1
}

actual=$(
  docker run --rm \
    --user "$(id -u):$(id -g)" \
    --volume "${repo_root}:/repo:ro" \
    --entrypoint /repo/bin/star_alignment_counts.py \
    "$star_image" \
    /repo/tests/fixtures/star/Log.final.out
)

[[ $actual == $'4\t7' ]] || {
  printf 'ERROR: configured STAR container returned an unexpected parser result\n' >&2
  exit 1
}

printf 'Configured STAR container parser contract passed.\n'
