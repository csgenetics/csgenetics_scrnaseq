---
name: csgenetics-scrnaseq
description: >
  Development guide for the CS Genetics public scRNA-seq pipeline
  (csgenetics/csgenetics_scrnaseq). Covers branching strategy, repository
  structure, testing, and contribution conventions. Use this skill whenever
  you are working in this repository — read it before making changes.
  Triggers: "how do I contribute", "branching strategy", "repo structure",
  "where do I start", "development workflow", "how to run tests",
  "which branch do I use".
---

# csgenetics_scrnaseq Development Guide

## What this repo is

The CS Genetics public scRNA-seq analysis pipeline. It processes droplet-based
single-cell RNA-Seq data from CS Genetics chemistry through QC, alignment, gene
annotation, cell calling, count matrix generation, and interactive HTML report
generation.

- **Repo:** `https://github.com/csgenetics/csgenetics_scrnaseq`
- **License:** MIT
- **Pipeline entry point:** `main.nf`
- **Pipeline config:** `nextflow.config`

## Stack

| Layer | Technology | Where to look |
|-------|-----------|---------------|
| Workflow orchestration | Nextflow (DSL2) | `main.nf`, `modules/` |
| Process scripts | Python | `bin/`, `modules/*/resources/usr/bin/` |
| Container images | Docker | `images/` (one subdirectory per tool image) |
| Containers config | Nextflow | `conf/images.config` |
| Execution platforms | See `nextflow.config` profiles | `conf/` |
| CI | CircleCI | `.circleci/config.yml` |
| Dependency management | Pixi | `pixi.toml`, `pixi.lock` |

## Repository structure

```
main.nf                  # Pipeline entry point
nextflow.config          # Global config, params, profiles, manifest (version)
conf/                    # Profile configs (aws, base, images, genomes, etc.)
modules/local/<name>/    # One directory per process (nf-core local-module layout)
bin/                     # Shared scripts (Python, gawk, plus compiled helpers)
tools/                   # Rust sources for the compiled helpers in bin/
templates/               # Jinja2 report template + Nextflow process templates
assets/vendor/           # Vendored JS/CSS/fonts, inlined so the report works offline
tests/nf/                # nf-test: whole-pipeline -stub DAG test
tests/python/            # pytest: report smoke test, number formatting
tests/regression/        # Output-equivalence comparator (strict + envelope modes)
images/                  # Docker image definitions (one dir per tool)
input_csv/               # Sample sheet examples
docs/                    # Documentation assets
.circleci/               # CircleCI workflow config
.claude/skills/          # Development skills (this file + others)
CHANGELOG.md             # Release notes -- update this for user-visible changes
```

## Branching strategy

| Branch | Purpose | Stability |
|--------|---------|-----------|
| `main` | Release branch | Most stable — only receives merges from `devel` |
| `devel` | Integration branch | Day-to-day development target |
| Feature branches | Individual changes | Created off `devel`, PR back to `devel` |
| Epic branches | Multi-ticket features | Created off `devel`, feature branches merge into epic, epic PRs to `devel` |

### Feature workflow

1. Create a feature branch off `devel`.
2. Implement the change, committing at logical checkpoints.
3. Open a PR targeting `devel`.
4. After review and CI, the PR is merged into `devel`.
5. `main` receives periodic release promotions from `devel`.

### Naming convention

- Feature branches: `GRT-XXX-short-descriptor` (referencing the tracking ticket)
- Epic branches: `epic/GRT-XXX-short-descriptor`

## Testing

CircleCI (`.circleci/config.yml`) runs four jobs on every PR:

| Job | What it does | Cost |
|-----|--------------|------|
| `check-no-internal-files` | Fails if agent working files, credentials or beast-local paths get committed. **This repo is public** — see `CLAUDE.md`. | seconds |
| `nf-test` | Whole-pipeline `-stub` DAG test plus per-process module tests. Catches wiring and output-declaration regressions without touching real data. | ~minutes |
| `report-smoke` | Renders the consolidated report in headless Chromium and asserts it is *functional*: no JS errors, plots drew SVGs, dropdown works, print reveals all panes. | ~minutes |
| `run-current-branch` | Launches a real end-to-end run on Seqera Platform against the `test` profile. | a full pipeline run |

Run the fast checks locally before pushing:

```bash
tests/check_no_internal_files.sh
NXF_VER=26.04.1 nf-test test tests/nf/pipeline_stub.nf.test
python -m pytest tests/python/ -v
nextflow run main.nf -profile test,docker      # full local run against test data
```

Available test profiles are defined in `nextflow.config` and `conf/test.config`. Note that the
`test` profile reads a **remote** sample sheet from `s3://csgx.public.readonly`; the copy at
`input_csv/test_input.csv` is for reference only.

For changes that could affect pipeline outputs, gate them with
`tests/regression/compare_outputs.py`, which has a strict (byte-identical) mode and an envelope
mode for the pipeline's inherent multimapper-ambiguity tolerance.

## Contribution guidelines

- All PRs target `devel`, never `main` directly.
- Follow existing code patterns and naming conventions.
- If adding a new Nextflow process, give it its own `modules/local/<name>/main.nf` and include a
  `stub:` block — the whole-pipeline `-stub` test depends on every process having one.
- If the process needs a new Docker image, add a Dockerfile under `images/` and pin the tag in
  `conf/images.config`. Images are pinned per code-version so older revisions keep working.
- Run the pipeline with a test profile before submitting a PR to verify nothing is broken.
- Record user-visible changes in `CHANGELOG.md`; bump `manifest.version` in `nextflow.config` for
  a release.
- **This repository is public.** Read the checklist at the top of `CLAUDE.md` before committing.

## Key files to read first

1. This skill file (you are here)
2. `CLAUDE.md` at repo root (if present) — coding standards and conventions
3. `main.nf` — pipeline structure
4. `nextflow.config` — parameters, profiles, global settings
5. `conf/images.config` — which Docker images are used by which processes
