# External scrnaseq NF26 Overhaul — Master Plan

Branch: `epic/external-nf26-overhaul` (off `devel`). Final PR target: `devel` (convention; merged by human).
Exploration findings: `plans/external_nf26_overhaul/exploration/01..06_*.md`.
Adversarial plan review (5 blockers/6 majors, all folded below): `plans/external_nf26_overhaul/plan_review.md`.

## Goal & invariants
Bring the customer-facing pipeline up to scratch on the latest Nextflow (26.04.1 installed), port
safe efficiency wins from the internal pipeline, build a robust test suite (nf-test + pytest + CI),
and improve report appearance — **without changing what the pipeline computes or its published outputs.**

Hard invariants (every phase re-verifies):
- Output-equivalence: the SET of published files, the `.metrics.csv` / `multisample_out.csv` schemas,
  and all metric VALUES are identical to the pre-change baseline. Baseline captured in Phase 0.
- Company name stays "CS Genetics" (Themis rename is a separate later PR).
- `main` and `devel` keep running: containers are pinned per code-version, so only the new code pulls
  any new image tags.

## Execution model (Seqera, not beast) — decided 2026-06-11
End-to-end / baseline / regression runs execute on **Seqera Platform** (the production path, what CI uses),
NOT locally on beast. Split:
- **Fast inner loop on beast:** `nf-test` (per-process) + `pytest` (`bin/` scripts) + the comparator harness
  (pulls run results from S3). No git push needed.
- **Seqera tier:** baseline, e2e, and regression runs in the **external workspace `40027684065767`** on the
  compute env **`ed0zjbIeIvKoUPuGsy8fA`** (`spot-us-east-1-scrnaseq-nf26_04-fusion`, Forge/aws-batch/us-east-1,
  Fusion v2 + Wave + NVMe, SPOT, created by mirroring internal primary `37Kx7gyrWg9NWK7X8lFBGy`). API base
  `https://api.cloud.seqera.io`, token `~/.auth_tokens .NextflowTower` (`pipelineuser`).
- **NF version** is pinned per launch (NOT per CE): baseline launches with `NXF_VER=25.10.2` (current code),
  candidate launches with `NXF_VER=26.04.1` — same mechanism the internal migration used.
- **Seqera launches from GitHub** (repo URL + revision), so running candidate code requires **pushing the
  epic/feature branch to `csgenetics/csgenetics_scrnaseq`** (user-authorised). Baseline runs current `devel`.
- **workDir** `s3://csg-tower-bucket`; **outdir** `s3://csg-nextflow/external_nf26_overhaul/<run>/`.

## Test data (reuse internal reduced sets — closes M6)
Rewrap internal's reduced CI FASTQs (`s3://csg-reference/circleci-reference-fastqs_v3/` human/mixed/mouse
~100k reads) into external 3-column `sample,fastq_1,fastq_2` CSVs; pair with external's reduced human chr1
index and the existing **reduced mixed reference** `s3://csgx.public.readonly/.../circleci-reference-fastqs_v3/mixed_reduced`.
This yields fast human + mixed-species + mouse benchmark/test profiles (the reduced mixed profile closes
review finding M6) runnable per-PR on Seqera. CE creds read `csg-reference` directly (no download process needed).

## Decisions (LOCKED 2026-06-11; refined post-review)
1. Module structure: **SPLIT** `modules/processes.nf` into per-process module files. Move only
   cleanly-isolated single-consumer scripts under module `resources/usr/bin/` (`moduleBinaries`); KEEP
   shared, `file()`-consumed, or cross-importing scripts in top-level `bin/` (`assign_multi_mappers.gawk`,
   the `create_multi_sample_report.py`↔`create_single_sample_report`/`cell_caller` import cluster). Inventory
   every `${baseDir}`/`${projectDir}` asset ref + intra-`bin/` import BEFORE moving anything. (B4, B5)
2. Report depth: **Option C — CONSOLIDATE** to one styled Jinja2 report (NOT the internal JS engine),
   removing/renaming the legacy per-sample + multi-sample HTMLs. **M2 resolved 2026-06-11: this is a BLESSED
   output-set change.** Metric VALUES stay identical (freeze `summary_statistics.py`); only presentation +
   the published-file SET change. Handle like the dedup re-baseline: enumerate removed/renamed/added files,
   re-cut the report published-file contract, and add a CHANGELOG + README customer deprecation note mapping
   old→new filenames.
3. Optimizations: **Tier-1 safe wins + Tier-2 (audited) + Rust fusion (GATED)**. STAR `--runThreadN` is
   **R1 not R0** (B2/B3): BAM record order feeds order-sensitive `gawk`/`umi_tools`, so it ships only after the
   order-perturbation test. The `umi_tools dedup` -> `count_table_builder` swap WILL change outputs; surface a
   quantified diff and get sign-off (deliberate re-baseline) before committing it.
4. Final PR target: **`devel`** (repo convention; epic branch -> devel, merged by the human, not the agent).

## Confirmed facts
- Current code FAILS to parse under NF26: `nextflow inspect main.nf -profile test` →
  `No such variable: sample_id` at `modules/processes.nf:498` (eager directive evaluation). So the Phase 0
  BASELINE must run under **Nextflow 25.10.2**, not NF26.
- External is already past most NF26 syntactic gates (env() quoted, lowercase channel, explicit
  closures, all `script:` no `shell:`, no `stub:`, no check_max, resourceLimits present).
- External has ZERO tests today; internal has a mature nf-test + pytest + CI suite to mirror.

## Phases

| # | Phase | Depends on | Credential |
|---|-------|-----------|-----------|
| 0a | **DONE 2026-06-11** — created external CE `ed0zjbIeIvKoUPuGsy8fA` (mirror internal primary, external `aws-tower` cred) | — | Seqera (have) |
| 0 | Baseline (regression §1) **ON SEQERA**: tag `baseline-pre-nf26`; launch `test` + reduced-mixed + 2 degenerate fixtures on CURRENT `devel` with `NXF_VER=25.10.2`; run twice (determinism probe); pull S3 results; archive `manifest.sha256` + `trace` | 0a | Seqera+AWS (have) |
| 1 | NF26 syntactic migration: publishDir/pattern closure-wrap (1 site inspect-confirmed @498, 15 by static audit — wrap + re-`inspect` iteratively until clean), remove `strict` flag, version flip, CI `NXF_VER`. Outputs: metrics EXACT, matrices within the **pre-determinism envelope** (see `baseline_determinism_finding.md`) | 0 | — |
| 1d | **Determinism: Path A (decided 2026-06-11, after investigation).** The pipeline is non-deterministic only in an irreducible multimapper ambiguity (~2 of 51k count entries, net-preserving, customer metrics byte-stable; STAR-seed forcing corrupted RSeQC, so rejected). NO code change — characterize the micro-envelope and gate future changes via `compare_outputs.py --envelope-max-flips N`: **metrics byte-exact + count-matrix per-barcode column sums byte-exact**. See `baseline_determinism_finding.md` | 1 | — |
| 2 | Module restructure (Decision 1): split `processes.nf` into per-process module files; shared scripts stay in `bin/`; inventory `${baseDir}` refs + intra-`bin/` imports first. Outputs == baseline | 1 | — |
| 3 | Tier-1 performance ONLY: `samtools -@`/`--write-index`, `featureCounts -T`, right-size cpus. **STAR `--runThreadN` gated behind the order-perturbation test (R1)**. Outputs == baseline; measure speedup | 2 | — |
| 3b | Tier-2 (consumer-audited R1) + gated Tier-3: Rust fusions (`bam-splitter`→`bam-assigner`) each behind an exact-equivalence nf-test; `umi_tools`→`count_table_builder` dedup swap = quantified + user-blessed re-baseline (regression §8). Sequenced AFTER Phase 5 freeze; re-cuts only dedup-downstream golden | 4, 5, sign-off | quay (have) |
| 4 | Test infra: `nf-test.config` + `tests/nf/` (must NOT includeConfig `conf/aws.config`), add `stub:` to every process, nf-tests in waves (`multi_sample_report` tested in the split wave), hashed S3 fixtures, pytest for `bin/`, CI guard scripts | 2,3 | AWS (have) |
| 5 | Report consolidation (Decision 2, Full Option C): build one consolidated styled Jinja2 report, remove/rename legacy HTMLs; freeze `summary_statistics.py` + golden metrics test (VALUES unchanged); re-cut the report published-file contract + CHANGELOG/README deprecation note (old→new filenames) | 1 | — |
| 6 | CI overhaul: add unit-tests + nf-tests jobs (`machine`, `resource_class: large`, `parallelism`; AWS creds as CircleCI env vars OR unsigned public-S3 fixtures); harden Tower e2e (pin `NXF_VER`, full-SHA `revision`, poll, validate S3 outputs) | 4 | AWS creds in CI (provision) |
| 7 | Seqera CI wiring (CE already built in 0a): create/refresh external `test` action bound to CE `ed0zjbIeIvKoUPuGsy8fA`; update CircleCI `TOWER_COMPUTE_ENV_ID`=`ed0zjbIeIvKoUPuGsy8fA` / confirm `TOWER_WORKSPACE_ID`=`40027684065767` out-of-band; apply full-SHA `revision` + `NXF_VER` preRunScript fixes; e2e run | 6 | Seqera (have) |
| 8 | Docker images: rebuild/push only changed images to `quay.io/csgenetics/*` (Tier-1 needs none); **verify quay push scope first**; pin tags | 3b,5 | quay push (verify) |
| 9 | Validation: local `test` (then heavier profiles), report preview, output-equivalence sign-off | all | AWS (have) |
| 9b | Docs/release: README runtime bump, CHANGELOG breaking-floor entry + `1.1.17` major-vs-minor decision, documented rollback (GAP-3 version flip = the revert commit) | 9 | — |
| 10 | Final review (3-agent adversarial panel) + open PR into `devel` (unmerged) | 9b | GitHub bot (have) |

## Parallelism (corrected per M1)
- Phase 1 first (the parse fix; unblocks everything).
- **Phases 2 and 3 are SERIAL, not parallel** — both mutate the process files; split (2) lands before
  threading (3). m3: also land the split before any process-editing so the CI test-mapping guard can work.
- Phase 5 (report) may run parallel to 2/3 once M2 is decided (it touches templates + report scripts, not
  the BAM/count path) — but it must freeze `summary_statistics.py` so Phase 3b's dedup re-baseline lands after.
- Phase 3b is gated: after Phase 4 (nf-test equivalence gates) + Phase 5 freeze + user sign-off.
- Phase 7 (Seqera) is credential/role-blocked and needs external-workspace provisioning; all else proceeds without it.

## Output-equivalence method (the linchpin — see regression_strategy.md, amended by review)
Phase 0 runs CURRENT code on the covering profiles under **NF 25.10.2** and archives per-file fingerprints by
output class (deterministic text / BAM / h5ad / plot-HTML), normalizing away volatile bytes. Key review fixes:
`bcGeneSummary.txt` compared as a SORTED (io_sequence, gene_name) multiset (not raw line order); `dedup.log`
normalized to the two integers `summary_statistics.py` consumes; a determinism probe that ALSO perturbs record
order (not just fixed-thread reruns) to expose any order-sensitive count drift before trusting an R0 claim.
Every later phase re-runs the covering profiles and diffs against this archive; any drift halts and is surfaced.

**Equivalence bar (Path A):** customer metrics CSVs, RSeQC, dedup counts, and the per-barcode count-matrix
column sums must be **byte-exact**; the count matrix's per-gene-per-barcode entries may differ only within the
measured net-preserving multimapper-ambiguity envelope (`--envelope-max-flips N`). A change to metrics or to any
column sum is a real regression. This holds high no-regression confidence without re-baselining any deliverable.
