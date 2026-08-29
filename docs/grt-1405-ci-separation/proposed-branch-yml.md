# Proposed GitHub check: PRs to `main` come only from `devel`

Paste-ready for Ben. The scrnaseq GitHub App cannot create workflow files
(`workflows` permission is deliberately absent). This file is the operator
apply for GRT-1405 deliverable 4.

Create `.github/workflows/branch.yml` on `devel` with the YAML below, merge
that to `main` so the workflow exists on the base, then mark the check
required on `main` branch protection.

Do not grant any agent `workflows` scope to land this.

## What it asserts

A pull request whose base is `main` is allowed only when:

1. the head repository is `csgenetics/csgenetics_scrnaseq` (not a fork), and
2. the head ref is exactly `devel`.

Anything else fails the check. A PR `feature -> main` fails. A PR from a
fork to `main` fails. A PR `devel -> main` passes.

The script is the whole job: if the `{ ... }` test is false, the step exits
non-zero and GitHub reports the check failed. There is no checkout and no
token use.

## Check name to mark required

**`PRs to main must come from devel`**

That string is both the workflow `name:` and the job `name:`. In the branch
protection dropdown it should appear as that. If GitHub shows a slashed
form, it is `PRs to main must come from devel / PRs to main must come from devel`.
The job id is `main-from-devel`.

The check will not appear in the dropdown until it has run once against
`main`. After the file is on `main`, open (or re-push) a PR targeting `main`
so the check fires, then select it.

## `pull_request_target`, not `pull_request`

Use `pull_request_target`. This is a public repository.

- `pull_request` runs the workflow file from the PR head. A PR that targets
  `main` from the wrong branch can include a rewritten `branch.yml` that
  keeps the same check name and always passes. On a fork PR, the fork's
  workflow file is what runs. The required check then gates nothing.
- `pull_request_target` always uses the workflow file from the **base**
  (`main`). The PR cannot rewrite or disable it.

`pull_request_target` defaults to a write-capable `GITHUB_TOKEN` and runs
in the base repository, which is how fork PRs have been turned into write
access. This workflow does not accept that default:

- `permissions: {}` grants nothing, including contents write
- there is no `actions/checkout` and no `run` that fetches PR code
- it only reads `github.event.pull_request.head` from the event payload

nf-core's `branch.yml` uses `pull_request`. That is the wrong event for a
**required** check on a public repo, because the PR can change the check.
We are not copying that part.

## YAML to paste as `.github/workflows/branch.yml`

```yaml
name: PRs to main must come from devel

on:
  pull_request_target:
    branches:
      - main

permissions: {}

jobs:
  main-from-devel:
    name: PRs to main must come from devel
    runs-on: ubuntu-latest
    steps:
      - name: Check PR source
        env:
          HEAD_REPO: ${{ github.event.pull_request.head.repo.full_name }}
          HEAD_REF: ${{ github.event.pull_request.head.ref }}
        run: |
          { [[ "$HEAD_REPO" == csgenetics/csgenetics_scrnaseq ]] && \
            [[ "$HEAD_REF" == devel ]]; }
```
