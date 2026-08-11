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
nextflow.config          # Global config, params, profiles
conf/                    # Profile configs (aws, base, images, genomes, etc.)
modules/                 # Nextflow DSL2 modules (process definitions)
bin/                     # Shared Python scripts
images/                  # Docker image definitions (one dir per tool)
templates/               # Nextflow process templates
input_csv/               # Sample sheet examples
docs/                    # Documentation assets
.circleci/               # CircleCI workflow config
.claude/skills/          # Development skills (this file + others)
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

The pipeline uses CircleCI for continuous integration. The CI config is at
`.circleci/config.yml` and runs end-to-end pipeline tests via Nextflow Tower
(Seqera Platform) against test data.

To run the pipeline locally with test data:

```bash
nextflow run main.nf -profile test,docker
```

Available test profiles are defined in `nextflow.config` and `conf/test.config`.

## Contribution guidelines

- All PRs target `devel`, never `main` directly.
- Follow existing code patterns and naming conventions.
- If adding a new Nextflow process, place it in the appropriate `modules/` subdirectory.
- If the process needs a new Docker image, add a Dockerfile under `images/`.
- Run the pipeline with a test profile before submitting a PR to verify nothing is broken.

## Key files to read first

1. This skill file (you are here)
2. `CLAUDE.md` at repo root (if present) — coding standards and conventions
3. `main.nf` — pipeline structure
4. `nextflow.config` — parameters, profiles, global settings
5. `conf/images.config` — which Docker images are used by which processes
