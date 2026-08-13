#!/usr/bin/env bash

# The source and its committed production binary must never drift. Build with
# the same reviewed Rust release recorded in the binary, run every Rust unit
# test, strip the static musl artifact, and compare the executable bytes.

set -euo pipefail

repo_root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
manifest="${repo_root}/tools/io_count_extract/Cargo.toml"
committed_binary="${repo_root}/bin/io_count_extract"
target_dir=$(mktemp -d)

cleanup() {
  if [[ -d $target_dir ]]; then
    find "$target_dir" -mindepth 1 -delete
    rmdir "$target_dir"
  fi
}
trap cleanup EXIT

[[ $(rustc --version) == 'rustc 1.93.0 (254b59607 2026-01-19)' ]] || {
  printf 'ERROR: io_count_extract must be verified with rustc 1.93.0\n' >&2
  exit 1
}

CARGO_TARGET_DIR="$target_dir" cargo test --locked --manifest-path "$manifest"
rustup target add x86_64-unknown-linux-musl
CARGO_TARGET_DIR="$target_dir" \
  cargo build --locked --release --target x86_64-unknown-linux-musl \
    --manifest-path "$manifest"

built_binary="${target_dir}/x86_64-unknown-linux-musl/release/io_count_extract"
strip "$built_binary"

if ! cmp --silent "$built_binary" "$committed_binary"; then
  printf 'ERROR: committed bin/io_count_extract does not match its reviewed source\n' >&2
  sha256sum "$built_binary" "$committed_binary" >&2
  exit 1
fi

printf 'io_count_extract source tests and committed binary parity passed.\n'
