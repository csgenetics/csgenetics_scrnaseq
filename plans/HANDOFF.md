# HANDOFF — external (customer-facing) scRNA-seq pipeline overhaul

> **Resuming after a disconnect?** Read this file first, then `git log --oneline -25`
> and `gh pr checks 78`. This is the single entry point; deeper detail is in the
> linked docs. Last updated 2026-06-18. Branch tip at write time: `abf2e54`.

---

## 1. One-paragraph status

We are overhauling the **customer-facing** scRNA-seq pipeline (`csgenetics/csgenetics_scrnaseq`)
on branch **`epic/external-nf26-overhaul`** → **PR #78 into `devel`** (open, **fully CI-green,
mergeable, intentionally unmerged** — it is left for human review/merge). The work has three
completed thrusts: (a) **Nextflow 26 migration** + per-process module split + nf-test harness;
(b) a **speedup campaign** (≈39% less compute, ~1.7x wall-clock on human, ~3.3x on mixed,
validated on real data with near-identical outputs); (c) a **consolidated HTML report** rebuild
plus a GUI-improvement pass. Everything is verified on real Seqera runs and/or in the exact CI
image. The only things outstanding are **user decisions** (see §4) and an optional deferred
dataviz polish (see §5). Nothing is in-flight; no Seqera runs are running.

## 2. The task & the hard invariants

- Bring the customer pipeline up to date on the latest Nextflow (26.04.1), port safe speedups
  from the internal pipeline, build a robust test suite, and improve the report — **without
  changing what the pipeline computes or its published outputs**, except where a change is
  explicitly *blessed* by the user.
- Company name stays **"CS Genetics"** (the Themis Bios rename is a separate later PR — do not
  do it here).
- `main`/`devel` must keep running: containers are pinned per code-version, so only the new code
  pulls new image tags.
- **Output-equivalence is the linchpin.** The pipeline has a tiny inherent multimapper-ambiguity
  envelope (a few of ~50k count-matrix entries can flip between two genes; per-barcode totals,
  cell calls, and summary metrics are stable). The user **accepted "very small changes of a count
  or so"** as the bar — NOT byte-identical. Gate every change with `tests/regression/compare_outputs.py`
  (strict + envelope modes), and on real data check `num_cells` (must hold) and read counts (within
  ~0.003%).

## 3. What's done (all on `epic/external-nf26-overhaul`, all verified)

**NF26 migration & structure**
- Parse-fix (closure-wrapped `publishDir`/`pattern`), `nextflowVersion = '!>=26.04.0'`, version `2.0.0`.
- 39 processes split from one `modules/processes.nf` into per-process `modules/local/<name>/` files.
  Shared scripts stay in top-level `bin/` (NOT moduleBinaries — `file()`/cross-imports would break).
- nf-test harness: `nf-test.config`, `tests/nf/`, `stub:` on every process, a whole-pipeline
  `-stub` DAG test, `tests/python/` pytest suite, and the regression comparator.

**Speedups (real-data validated; full detail in `plans/external_speedup/CAMPAIGN_LOG.md`)**
- `io_count`: BusyBox-awk → static Rust binary `bin/io_count_extract` (~150x; byte-identical).
- `dedup`: split-by-contig parallel `umi_tools` via `bin/dedup_by_contig.sh` (counts identical;
  only the intermediate `dedup.bam` representative reads differ).
- `multimapper_*` + filters + `initial_feature_count`: deterministic canonical assignment + threaded
  `samtools` (the deterministic-canonical change is a one-time ~0.3% count shift; cell calls unchanged).
- Resource right-sizing (STAR `runThreadN` pinned to original **8**; cpus right-sized for throughput).
- **Reverted dead-ends (measured, not guessed):** STAR thread-bumping (memory-bound, no gain);
  name-hash multimapper split (2.3x *slower* on heavy data); RSeQC contig-split (slower on real BAMs).
  Lesson: parallel ≠ faster — measure on real/heavy data before keeping.

**Bugs found by heavy + multi-lane testing (all fixed, pre-existing / also on `main`)**
- `filter_count_matrix` channel race on Seqera/Fusion (`order_integer_first` → replaced with `.join`).
- Multi-lane fan-out crashing `qc_cascade_plot_multi` (`.unique()` the per-lane threshold channel).
- Empty-sample dedup log wrote a literal `\n` (→ `printf`).

**Consolidated report + GUI pass**
- New `consolidated_report` process/`bin/create_consolidated_report.py` replacing the old single/multi
  report processes; offline-safe (Plotly/Bootstrap/font inlined); searchable sample dropdown.
- Fixed: blank plots (Plotly `<script>` was after the `newPlot()` calls → moved to `<head>`); dropdown
  collapsing to one option; mixed-species report (cell metrics classified per-metric, not "Cell metrics").
- GUI improvements landed: thousands-separator display formatting (CSV stays raw/byte-identical),
  removed the misleading small-N violin plot, print/PDF stylesheet (was silently dropping all-but-first
  sample), **Jinja autoescape hardening** (sample ids are customer input — XSS-safe now),
  **WCAG 2.1/2.2 AA accessibility** (zero axe violations), and **run-provenance** header/footer/key-results
  (threaded `workflow.*` + params via `main.nf` → process → generator → template; non-judgemental).
- **CI smoke test** (`tests/python/test_report_smoke.py` + `report-smoke` CircleCI job): renders a fixture
  report in headless Chromium and asserts it is *functional* (0 JS errors, plots drew SVGs, dropdown
  works, print reveals all panes, CSV separator-free, XSS-escaped). This guards the exact class of bug
  that shipped twice.
- Deleted dead `bin/create_multi_sample_report.py` (no process used it after consolidation).

**CI:** `nf-test`, `report-smoke`, `run-current-branch`, `claude-review` all **pass** on the tip.

## 4. Open items — these are the USER's decisions (not blockers)

1. **Review / merge PR #78.** It is green and mergeable; `mergeStateStatus=BLOCKED` is only the
   required human-approval gate. The agent does not merge it.
2. **Multimapper memory headroom:** on heavy samples `multimapper_assignment` OOMs at the 12GB base
   reservation and self-heals via Nextflow retry at 24GB (≈20–28 min wasted on the doomed first
   attempt). Bumping the base 12→24GB trades that wall-clock back for worse instance packing on the
   common (light) case. Flagged, **not acted on** — it's a cost/packing tradeoff for the user.
3. **Docker images:** the report uses the existing `quay.io/csgenetics/html_build:0.1.0` (no new image
   needed). The Rust binaries are committed static musl binaries (run in existing public containers).
   If a future change needs a new/updated quay image, that's a push the user must authorise.

## 5. Deferred (well-specified, not done)

- Report dataviz improvements flagged by the 8-agent review: cell-caller **log** count axis,
  QC-cascade **% retained** labels, a **barcode-rank knee** plot. (Some need new computation.)
- Optional commercial **methods / tool-versions appendix** (from pinned `conf/images.config` tags).
- The user explicitly **declined**: verdict/pass-fail reporting (scientist's domain), cross-sample
  **column** sorting (metrics are in rows — meaningless), Plotly **minification** (future features may
  need removed traces), **mobile** layout work. Support email + marketing URL omitted (don't exist yet).

## 6. How to resume specific kinds of work

- **Run / re-run on Seqera** (this is the production path; do NOT run locally on beast for validation):
  external workspace `40027684065767`, compute env `ed0zjbIeIvKoUPuGsy8fA`
  (`spot-us-east-1-scrnaseq-nf26_04-fusion`), API base `https://api.cloud.seqera.io`, token
  `~/.auth_tokens` key `NextflowTower`. Seqera launches from the **GitHub** repo + revision, so push
  the branch first. NF version is pinned **per launch** (`preRunScript: export NXF_VER=26.04.1`).
  Use a **full SHA** as `revision` (short SHA fails to fetch on non-default branches).
- **Report-only re-render via RESUME** (cheap): `POST /workflow/launch` with `resume:true` +
  `sessionId:<orig run's session>` + `revision:<new tip>` + same params/outdir → only
  `consolidated_report` re-runs, everything heavy is cached. (Used for the report fixes; see RESUME memory.)
- **Verify a report — ALWAYS render it, never grep:** headless Chromium screenshot +
  `--enable-logging=stderr` for JS console errors, and Playwright (installed) for interactions/keyboard.
  axe-core at `/nssd2/humebc/webapp/cs-genetics-webapp/node_modules/axe-core/axe.min.js`. The CI smoke
  test (`tests/python/test_report_smoke.py`) is the gate; reproduce CI in Docker with `cimg/python:3.12`.
- **Output-equivalence:** `tests/regression/compare_outputs.py` (strict + `--envelope-max-flips N`).
- **Real validation data:** internal datasets are at `/nssd2/humebc/pipeline-runs/fastqs/`
  (KOL0054 mixed, MOR034/MOR036 human); heavy lane-merged copies + last validation outputs at
  `s3://csg-nextflow/external_nf26_overhaul/validation/`. Long-term-bucket FASTQs are mostly Glacier.

## 7. Credentials / infra cheat-sheet (all present on beast, account B)

| Need | Where |
|---|---|
| AWS S3 (read csg-reference etc.) | `AWS_SHARED_CREDENTIALS_FILE=/nssd2/humebc/agent-secrets/scrnaseq/aws/credentials`, `AWS_CONFIG_FILE=.../aws/config` (user `scrnaseq-agent`) |
| Seqera Platform token | `~/.auth_tokens` → `.NextflowTower` (user `pipelineuser`) |
| GitHub bot (PRs/comments) | `GH_TOKEN=$(cat /nssd2/humebc/agent-secrets/scrnaseq/github/token)` (bot `csgenetics-scrnaseq-agent`) |
| CircleCI API (read CI logs — user-authorised 2026-06-18) | `~/.auth_tokens` → `.CircleCI` |
| quay.io | docker already logged in (pull works; push scope unverified) |
| Token usage (5h/7d) | sidecar socket `/nssd2/humebc/agent-secrets/claude-usage-sidecar/sub-b/claude-usage.sock` → `http://x/usage` (account B) |
| Hardware | beast: 128 CPU / 251 GB RAM; Docker 27, Singularity 3.8.5, nextflow + nf-test on PATH |

## 8. Deeper docs & memory (read on demand)

- `plans/external_speedup/CAMPAIGN_LOG.md` — full speedup campaign log (every strategy, measured numbers, reverts).
- `plans/external_nf26_overhaul/` — master_plan, plan_review (adversarial), regression_strategy,
  baseline_determinism_finding, optimization_assessment, report_consolidation.
- Memory notes: `/nssd2/humebc/.claude-b/projects/-nssd2-humebc-external-csgenetics-scrnaseq/memory/`
  (`RESUME-current-state.md` is the rolling live state; `MEMORY.md` is the index).
- Repo rules: `CLAUDE.md` (fail loud, no defensive programming, check file content in reviews, verify
  claims with line numbers, `.nextflow.log`-based work-dir lookup).

## 9. Hard-won lessons (don't relearn these)

- **Render reports to verify them.** Grepping HTML for tokens is NOT verification — two report bugs
  shipped that way. Screenshot + console + Playwright.
- **Reproduce CI in Docker before pushing a CI fix** (`cimg/python:3.12`) — saves red-CI ping-pong.
- **Optimise wall-clock (time-to-result), not just CPU-minutes** — they diverge; measure on heavy/real
  data. STAR/featureCounts threading gave nothing at our read counts; the real cost was BusyBox awk +
  single-threaded sorts.
- **The pipeline is non-deterministic at the count-matrix edge by nature** (multimapper ambiguity);
  `PYTHONHASHSEED=0` fixed the umi_tools part, the gawk canonical-representative change fixed another.
  Forcing full determinism via STAR flags corrupts RSeQC — don't.
- **Heavy + multi-lane testing finds what single-lane smoke tests miss** (the three pre-existing bugs).
