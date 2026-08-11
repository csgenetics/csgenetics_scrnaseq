#!/usr/bin/env bash

# Focused regression tests for the public-repository guard. Each case builds a fresh
# temporary Git index because the guard intentionally scans staged content.

set -euo pipefail

repo_root=$(git rev-parse --show-toplevel)
scratch_root=$(mktemp -d)
trap 'rm -rf -- "$scratch_root"' EXIT
success_message='OK: no internal-only paths, credentials, private S3 endpoints, local paths or external Conda sources staged.'

create_case_repo() {
  local name=$1
  local case_repo="$scratch_root/$name"

  mkdir -p "$case_repo/tests"
  cp "$repo_root/tests/check_no_internal_files.sh" "$case_repo/tests/"
  git -C "$case_repo" init -q
  git -C "$case_repo" add tests/check_no_internal_files.sh
  printf '%s\n' "$case_repo"
}

expect_pass() {
  local name=$1
  local content=$2
  local case_repo output

  case_repo=$(create_case_repo "$name")
  printf '%s\n' "$content" > "$case_repo/fixture.txt"
  git -C "$case_repo" add fixture.txt
  if ! output=$(cd "$case_repo" && tests/check_no_internal_files.sh 2>&1); then
    printf 'FAIL: expected %s to pass\n%s\n' "$name" "$output" >&2
    return 1
  fi
  if [[ "$output" != "$success_message" ]]; then
    printf 'FAIL: %s returned unexpected output\n%s\n' "$name" "$output" >&2
    return 1
  fi
}

expect_recipe_path_pass() {
  local name=$1
  local content=$2
  local case_repo output

  case_repo=$(create_case_repo "$name")
  mkdir -p "$case_repo/conda-recipes/example"
  printf '%s\n' "$content" > "$case_repo/conda-recipes/example/meta.yaml"
  git -C "$case_repo" add conda-recipes/example/meta.yaml
  if ! output=$(cd "$case_repo" && tests/check_no_internal_files.sh 2>&1); then
    printf 'FAIL: expected %s to pass\n%s\n' "$name" "$output" >&2
    return 1
  fi
  if [[ "$output" != "$success_message" ]]; then
    printf 'FAIL: %s returned unexpected output\n%s\n' "$name" "$output" >&2
    return 1
  fi
}

expect_recipe_path_failure() {
  local name=$1
  local content=$2
  local case_repo output

  case_repo=$(create_case_repo "$name")
  mkdir -p "$case_repo/conda-recipes/example"
  printf '%s\n' "$content" > "$case_repo/conda-recipes/example/meta.yaml"
  git -C "$case_repo" add conda-recipes/example/meta.yaml
  if output=$(cd "$case_repo" && tests/check_no_internal_files.sh 2>&1); then
    printf 'FAIL: expected %s to fail\n' "$name" >&2
    return 1
  fi
  if [[ "$output" != *"Conda recipe source path does not resolve inside this public repository"* ]]; then
    printf 'FAIL: %s failed without the expected diagnostic\n%s\n' "$name" "$output" >&2
    return 1
  fi
  if [[ "$output" == *"$content"* ]]; then
    printf 'FAIL: %s repeated the recipe content in its diagnostic\n' "$name" >&2
    return 1
  fi
}

expect_failure() {
  local name=$1
  local content=$2
  local diagnostic=$3
  local case_repo output

  case_repo=$(create_case_repo "$name")
  printf '%s\n' "$content" > "$case_repo/fixture.txt"
  git -C "$case_repo" add fixture.txt
  if output=$(cd "$case_repo" && tests/check_no_internal_files.sh 2>&1); then
    printf 'FAIL: expected %s to fail\n' "$name" >&2
    return 1
  fi
  if [[ "$output" != *"$diagnostic"* ]]; then
    printf 'FAIL: %s failed without the expected diagnostic\n%s\n' "$name" "$output" >&2
    return 1
  fi
  if [[ "$output" == *"$content"* ]]; then
    printf 'FAIL: %s repeated sensitive fixture content in its diagnostic\n' "$name" >&2
    return 1
  fi
}

expect_staged_failure() {
  local staged_content=$1
  local working_content=$2
  local case_repo output

  case_repo=$(create_case_repo staged_content_wins)
  printf '%s\n' "$staged_content" > "$case_repo/fixture.txt"
  git -C "$case_repo" add fixture.txt
  printf '%s\n' "$working_content" > "$case_repo/fixture.txt"

  if output=$(cd "$case_repo" && tests/check_no_internal_files.sh 2>&1); then
    printf 'FAIL: expected staged content to be scanned\n' >&2
    return 1
  fi
  if [[ "$output" != *"private S3 endpoint literal staged"* ]]; then
    printf 'FAIL: staged-content case lacked the expected diagnostic\n%s\n' "$output" >&2
    return 1
  fi
  if [[ "$output" == *"$staged_content"* ]]; then
    printf 'FAIL: staged-content case repeated the endpoint in its diagnostic\n' >&2
    return 1
  fi
}

scheme_separator=':/'
scheme_separator+='/'
s3_scheme="s3${scheme_separator}"
s3a_scheme="s3a${scheme_separator}"
s3n_scheme="s3n${scheme_separator}"
https_scheme="https${scheme_separator}"
uppercase_s3_scheme=${s3_scheme^^}

private_bucket='company-results'
public_bucket='research.public.readonly'
public_suffix_bucket='community-public'
example_bucket='example-bucket'
aws_s3_service='s3'
aws_domain='amazonaws.com'

expect_pass public_bucket \
  "${s3_scheme}${public_bucket}/catalog/genes.gtf"
expect_pass public_suffix \
  "${s3a_scheme}${public_suffix_bucket}/catalog/genes.gtf"
expect_pass example_bucket \
  "${s3n_scheme}${example_bucket}/documentation/file.fastq.gz"
expect_pass public_virtual_host \
  "${https_scheme}${public_bucket}.${aws_s3_service}.${aws_domain}/catalog/genes.gtf"
expect_pass example_path_style \
  "${https_scheme}${aws_s3_service}.eu-west-2.${aws_domain}/${example_bucket}/file.fastq.gz"
expect_pass unrelated_https \
  "${https_scheme}example.com/${aws_s3_service}/${private_bucket}/documentation"
expect_pass injected_private_location \
  'results_uri=${PIPELINE_OUTDIR_ROOT%/}/run'

expect_failure private_s3 \
  "${s3_scheme}${private_bucket}/customer/run.csv" \
  "private S3 endpoint literal staged"
expect_failure private_s3a \
  "${s3a_scheme}${private_bucket}/customer/run.csv" \
  "private S3 endpoint literal staged"
expect_failure private_s3n \
  "${s3n_scheme}${private_bucket}/customer/run.csv" \
  "private S3 endpoint literal staged"
expect_failure uppercase_scheme \
  "${uppercase_s3_scheme}${private_bucket}/customer/run.csv" \
  "private S3 endpoint literal staged"
expect_failure private_virtual_host \
  "${https_scheme}${private_bucket}.${aws_s3_service}.${aws_domain}/customer/run.csv" \
  "private S3 endpoint literal staged"
expect_failure private_regional_virtual_host \
  "${https_scheme}${private_bucket}.${aws_s3_service}.eu-west-2.${aws_domain}/run.csv" \
  "private S3 endpoint literal staged"
expect_failure private_legacy_regional_virtual_host \
  "${https_scheme}${private_bucket}.${aws_s3_service}-eu-west-2.${aws_domain}/run.csv" \
  "private S3 endpoint literal staged"
expect_failure private_china_virtual_host \
  "${https_scheme}${private_bucket}.${aws_s3_service}.cn-north-1.${aws_domain}.cn/run.csv" \
  "private S3 endpoint literal staged"
expect_failure private_path_style \
  "${https_scheme}${aws_s3_service}.${aws_domain}/${private_bucket}/customer/run.csv" \
  "private S3 endpoint literal staged"
expect_failure private_regional_path_style \
  "${https_scheme}${aws_s3_service}.eu-west-2.${aws_domain}/${private_bucket}/run.csv" \
  "private S3 endpoint literal staged"
expect_failure private_legacy_regional_path_style \
  "${https_scheme}${aws_s3_service}-eu-west-2.${aws_domain}/${private_bucket}/run.csv" \
  "private S3 endpoint literal staged"

expect_staged_failure \
  "${s3_scheme}${private_bucket}/staged/run.csv" \
  "${s3_scheme}${example_bucket}/working-tree/run.csv"

expect_recipe_path_pass in_repository_recipe_source \
  $'source:\n  path: ../../images/qc'
expect_recipe_path_pass quoted_in_repository_recipe_source \
  $'source:\n  - path: "../../images/qc" # tracked source'
expect_recipe_path_failure escaping_recipe_source \
  $'source:\n  path: ../../../sibling-project/images/qc'
expect_recipe_path_failure absolute_recipe_source \
  $'source:\n  path: /outside/repository/images/qc'
expect_recipe_path_failure dynamic_recipe_source \
  $'source:\n  path: {{ environ["QC_SOURCE"] }}'
expect_recipe_path_failure shell_variable_recipe_source \
  $'source:\n  path: $QC_SOURCE/images/qc'

access_key='AK'
access_key+='IA'
access_key+='ABCDEFGHIJKLMNOP'
expect_failure credential "$access_key" "possible credential staged"

printf 'OK: public-repository guard regression tests passed.\n'
