# 04 - Test Infrastructure Map (Internal reference -> External blueprint)

Scope: map INTERNAL (`/nssd2/humebc/internal/rnaseq`) test infrastructure and define a
concrete BLUEPRINT for building an equivalent suite in EXTERNAL
(`/nssd2/humebc/external/csgenetics_scrnaseq`). External currently has ZERO test
infrastructure.

---

## 1. Internal nf-test setup

### Config (`internal/rnaseq/nf-test.config`)
```
testsDir "modules"
workDir   System.getenv("NFT_WORKDIR") ?: ".nf-test"   # per-process workdir isolation for parallel runs
configFile "tests/nf/nextflow.config"
plugins { load "nft-utils@0.0.9" }
triggers "nextflow.config", "nf-test.config", "tests/nf/nextflow.config"
```
- `NFT_WORKDIR` indirection is load-bearing: lets the parallel triage loop give each test
  its own workdir to avoid `NoSuchFileException` from `.nf-test/` contention.

### Directory convention
- Process tests: `modules/<module>/tests/<process>.nf.test` (+ committed `<process>.nf.test.snap`).
- Pipeline test: `tests/nf/pipeline_stub.nf.test` (uses `nextflow_pipeline`).
- nf-test nextflow config: `tests/nf/nextflow.config`.
- Snapshot exclusion: `tests/nf/.nftignore`.
- Per-test fixtures committed in-repo: `tests/nf/fixtures/` (e.g. `human.100000.test.non_pbmc.input_csv.csv`).
- Scale: 210 `.nf.test` files across `modules/*/tests/`; 30 committed `.snap` files (snapshots
  used only where output is deterministic text). ~75 processes / 15 modules.

### nf-test nextflow.config (`tests/nf/nextflow.config`)
- `includeConfig '../../conf/images.config'` (real containers), `plugins { id 'nf-amazon' }`.
- Fixtures live on S3, addressed via params:
  - `params.nftest_data_base = "s3://csg-reference/internal_nf_tests_data"`
  - `params.nftest_references_base = "s3://csgx.public.readonly/resources/circleci-reference-fastqs_v3/mixed_reduced"`
- AWS profile selection is a top-level ternary (NF26 config parser v2 forbids `if`): honor
  `AWS_ACCESS_KEY_ID`/`AWS_PROFILE` env (CI), else `[csg-rnaseq-nftest]` named profile on beast.
- Docker enabled with `runOptions = '-u $(id -u):$(id -g)'`; `nextflow.enable.moduleBinaries = true`
  (required for script-based processes whose scripts live in `resources/usr/bin/`).
- Retry policy `errorStrategy 'retry'`, `maxRetries 2`; resources scale with `task.attempt`.
- A full block of pipeline params is pre-set so process tests do not need a real run context
  (species prefixes, purity, dedup, barcode pattern, mito/ribo strings, etc.).

### Fixture strategy
- Default fixtures are MIXED-species (human+mouse) to exercise the superset of code paths;
  single-species cases added only where logic genuinely diverges (`cell_caller`, `annotate_cells`).
- Only the `star` test downloads the STAR index (~200-400 MB); all downstream tests consume
  pre-aligned BAM fixtures -> ~90% of tests stay fast/independent.
- Fixtures referenced via string concatenation, never GString interpolation:
  `file(params.nftest_data_base + '/path', checkIfExists: true)` (nf-test evaluates `"${...}"`
  before NF runs, so interpolation resolves to null).
- CONTENT-HASH-IN-FILENAME is mandatory and CI-enforced: `sample1_R1.a5d7bf18.fastq.gz`.
  Default = 8-hex hash in filename before the semantic suffix. Exception: when a production
  script parses metadata from the filename (glob/split), hash goes in a parent DIRECTORY
  instead (`barcode_diversity/4a5cb186/sample1.1412380.barcode_frequencies.tsv`).

### Snapshot discipline
- Snapshot ONLY deterministic text: CSV metrics, JSON, `env()` counts. Snap stores md5
  checksums + filenames (e.g. count_reads snap: `[["sample1","100000","sample1_R1...:md5,..."]]`).
- NEVER snapshot binary (BAM, h5ad) or nondeterministic (HTML, PNG) outputs -> assert
  `.exists()` and `.size() > N` instead.
- Snapshots must be generated locally (`nf-test test <file>`) and COMMITTED; CI runs with
  `--ci` which refuses to auto-create missing snapshots.
- Intentional change: `nf-test test <file> --update-snapshot`, then review `git diff *.snap`.

### Stub tests
- Every process test pairs a positive case with a stub case (`options "-stub"`), and every
  process has a `stub:` block in `main.nf` producing all declared outputs (2-10 lines).
- Stub assertions: bare `assert process.success` OR `assertAll({ assert process.success })`.
- Stub blocks enable the pipeline-level DAG validation test (`tests/nf/pipeline_stub.nf.test`,
  `options "-stub -profile ci_test"`) asserting `workflow.success`,
  `workflow.trace.succeeded().size() > 0`, `workflow.trace.failed().size() == 0`.

### Tagging convention
- Every process test tags `"modules"` plus the module dir name (e.g. `"qc"`, `"gene_annotation"`).
- Pipeline test tags `"pipeline"`, `"stub"`. Tier tags (`tier:core`) are aspirational, not in use.

### Representative tests read
- `modules/qc/tests/count_reads.nf.test` — env-count output + snapshot + stub case.
- `modules/gene_annotation/tests/initial_feature_count.nf.test` — BAM output, existence/size
  asserts (no snapshot), stub case.
- `modules/qc_reporting/tests/qc_stats_multi_sample.nf.test` — multi-sample list input, CSV
  existence/size asserts, stub case.
- `tests/nf/pipeline_stub.nf.test` — pipeline DAG + a regression case (PBMC gating) driven
  purely by an in-repo input_csv fixture and `workflow.success`.

---

## 2. Internal pytest setup

### Config (`internal/rnaseq/pytest.ini`)
- `testpaths = tests`, `python_files = test_*.py`, `addopts = -v --strict-markers --tb=short`.
- Markers: `unit`, `integration`, `regression`.

### Layout (`internal/rnaseq/tests/`)
- Mirrors module structure: `pipeline_report/`, `metrics_collections/`, `internal_metrics/`,
  `stress/`, `barcode_diversity/`, `summary_reporting/`, `pbmc_reference_check/`, `utilities/`,
  `full/`, plus top-level `test_interp_extrap_module.py`.
- ~45 custom `test_*.py` files (the 4265 raw count is dominated by vendored `.pixi` env tests —
  ignore those). Heavy concentration in `pipeline_report/unit/` (~30 files) and
  `barcode_diversity/` (unit + integration incl. Rust-binary, cross-validation, end-to-end).
- Tests organized `unit/` vs `integration/` per module.
- `conftest.py` per module sets `sys.path` to production code (`bin/`, `bin/utilities/`,
  `resources/usr/bin/...`); root `conftest.py` only sets an anndata write flag. No manual
  PYTHONPATH needed. Shared assets in `tests/fixtures/`.

### Environments (pixi)
- Two pixi manifests under `tests/requirements/`:
  - `tests/requirements/pipeline_report/pixi.toml` — pipeline_report / metrics_collections /
    internal_metrics tests (lighter).
  - `tests/requirements/full/pixi.toml` — adds scanpy/matplotlib/pysam for stress +
    barcode_diversity; runs the whole suite.
- Must run via pixi (`cd tests/requirements/<env> && pixi run pytest tests/... -v`) to pin
  Python 3.11 + NumPy 2.x. pytest.ini at repo root drives discovery.

### Coverage character
- Pure-Python unit tests of the report/metrics/stress/barcode-diversity LOGIC (the contents of
  `bin/` and `modules/*/resources/usr/bin/` scripts), independent of Nextflow.
- Plus a JS suite (Vitest+jsdom) for the pipeline report UI: `modules/pipeline_report/js/`
  (37 files, 1637+ tests) — relevant only if external adopts a JS-based report.

---

## 3. CI design

### Internal (`internal/rnaseq/.circleci/config.yml`) — 3 workflows
1. `unit-tests` -> `run-unit-tests` (docker `cimg/python:3.11-node`): installs pixi, runs
   `pytest tests` via `tests/requirements/full/pixi.toml`, plus JS tests
   (`npm ci`/`npm test` in `modules/pipeline_report/js`).
2. `nf-tests` -> `run-nf-tests` (machine executor, `resource_class: large`, **parallelism: 4**,
   `NXF_VER=26.04.1`): installs Java17 + Nextflow + nf-test; on shard 0 runs
   `.circleci/check_nf_tests.sh` (modified `modules/*/main.nf` must have tests; `|| true` = warn)
   and `.circleci/check_nf_fixture_hashes.sh` (HARD fail; enforces fixture content hashes +
   a registry of filename-parsing scripts); then
   `nf-test test --ci --verbose --shard $((IDX+1))/$TOTAL --junitxml nf-test-results.xml`
   with `store_test_results`.
3. `e2e-tests` -> `e2e-mixed-species` + `e2e-single-species`: launch on Seqera Tower via REST
   (create action -> launch with `species_profile` test params -> poll status -> validate
   expected S3 output files: `.pipeline_report.html`, `metrics_definitions.html`,
   `multiqc/multisample_multiqc.html`, `extrapolation_interpolation/interp_extrap.predictions.csv`
   -> delete action). Reusable commands: `tower-create-action`, `tower-poll-status`,
   `tower-validate-s3-outputs`.

### CI helper scripts
- `check_nf_tests.sh`: `git diff --name-only origin/main...HEAD -- 'modules/*/main.nf'`, every
  modified module must have a `tests/*.nf.test`. (Internal layout: one module per dir.)
- `check_nf_fixture_hashes.sh`: greps `params.nftest_data_base + '...'` out of every test,
  requires an 8-hex hash in filename or parent dir; auto-discovers glob/split filename parsers
  in scripts and requires them registered in `GLOB_SENSITIVE_SUFFIXES` (with a
  `REVIEWED_SAFE_SCRIPTS` allowlist). This is the most sophisticated guardrail.

### External (`external/.circleci/config.yml`) — TODAY
- ONE workflow `nextflow-tower` / one job `run-current-branch`: create Tower action -> launch
  with `configProfiles:["test"]` -> monitor status -> delete endpoint. NO unit tests, NO
  nf-tests, NO output validation, NO sharding, NO `NXF_VER` pin (uses Tower default).
- `.github/workflows/` only holds `claude.yml` + `claude-code-review.yml` (bot, not tests).

---

## 4. External current state (what exists to test)

- NO `nf-test.config`, NO `tests/`, NO `pytest.ini`, NO `*.nf.test`, NO `conftest.py`
  (the only matches are inside vendored `.pixi/.../nf_core` and site-packages — ignore).
- Processes: **39 processes in a single monolithic `modules/processes.nf`** (NOT one-dir-per-module
  like internal). Key processes (line in processes.nf): `merge_lanes`(131), `qc`(158),
  `star`(191), `gtf2bed`(249), `run_rseqc`(266), `initial_feature_count`(291),
  `filter_for_UMRs_mismatch`(319), `umr_transcript_assignment`(335), `umr_exon_assignment`(351),
  `filter_for_multimappers_mismatch`(377), `multimapper_*`(392/435), `merge_*_bams`(465/480/495),
  `count_high_conf_annotated_umr_multimap`(512), `single_sample_multiqc`(532),
  `multi_sample_multiqc`(559), `sort_index_bam`(586), `dedup`(607), `io_count`(642),
  `count_matrix`(674), `cell_caller`(705), `filter_count_matrix`(733), `categorize_reads`(764),
  `summary_statistics`(790), `qc_cascade_plot_single`(820), `qc_cascade_plot_multi`(843),
  `single_summary_report`(863), `multi_sample_report`(885), plus download/setup processes.
- ZERO `stub:` blocks in `processes.nf` (grep count 0) -> stub tests + pipeline DAG test are
  not yet possible without adding stubs.
- Python/awk/gawk logic lives in top-level `bin/` (categorize_reads.py, cell_caller.py,
  count_matrix.py, filter_count_matrix.py, features_names.py, summary_statistics.py,
  qc_cascade_plot.py, create_single/multi_sample_report.py, create_curated_barcode_correction_list.py,
  *.awk/*.gawk) — directly pytest-able.
- Test profiles already exist: `test`, `test_pbmc_4_sample_full`, `test_hsap_mmus_2_sample_full`,
  `test_singularity`, `test_conda` (nextflow.config:93-148), each layering
  `conf/base.config` + `conf/test*.config` + `conf/images.config`. Good seed for nf-test and
  pipeline-stub profiles.

---

## 5. BLUEPRINT for the external test suite

### A. nf-test scaffolding (port internal conventions)
1. Add `nf-test.config` at external root. Decision needed on `testsDir`: internal uses
   `testsDir "modules"`, but external has a single `modules/processes.nf`. Recommended:
   `testsDir "tests/nf/modules"` (tests grouped by process file name there) OR split
   `processes.nf` into per-module dirs during the overhaul and mirror internal exactly.
   Keep `workDir System.getenv("NFT_WORKDIR") ?: ".nf-test"` and `configFile
   "tests/nf/nextflow.config"`.
2. Add `tests/nf/nextflow.config`: `includeConfig '../../conf/images.config'`,
   docker enabled with uid/gid runOptions, retry policy, `nextflow.enable.moduleBinaries`
   only if external moves script staging into `resources/usr/bin/` (today scripts are in `bin/`
   — see open questions). Define `params.nftest_data_base` / `params.nftest_references_base`
   pointing at an external S3 fixture prefix (reuse `s3://csgx.public.readonly/...mixed_reduced`
   references that already feature in the existing Tower CI input).
3. Add `tests/nf/.nftignore`.

### B. Which processes need nf-tests (priority waves, external names)
- Wave 1 (shell-only, prove framework): `merge_lanes`, `io_count`, `sort_index_bam`,
  `filter_for_UMRs_mismatch`, `count_high_conf_annotated_umr_multimap`.
- Wave 2 (binary/heavy): `qc` (Rust-in-Docker), `star` (needs index — only test that pulls it),
  `initial_feature_count` (featureCounts), `run_rseqc`/`gtf2bed`.
- Wave 3 (Python script processes): `count_matrix`, `cell_caller`, `filter_count_matrix`,
  `categorize_reads`, `summary_statistics`.
- Wave 4 (reporting): `single_sample_multiqc`, `multi_sample_multiqc`,
  `qc_cascade_plot_single/multi`, `single_summary_report`, `multi_sample_report`.
- Wave 5: remaining merge/annotation/dedup/download processes.
- Each test: positive case + (param-dependent) alt case + stub case. ADD a `stub:` block to
  every process in `processes.nf` first (currently 0) so stub tests + the pipeline DAG test work.
- Snapshot policy identical to internal: snapshot deterministic CSV/JSON/env counts
  (`io_count`, `summary_statistics` metrics, `count_matrix` summary CSVs); assert existence/size
  for BAM/h5ad/HTML/PNG.

### C. Fixtures / test data
- Stand up an external S3 fixture prefix mirroring `internal_nf_tests_data/` layout
  (`fastq/mixed/`, `process_fixtures/<group>/`, `references/`). Reuse the existing public
  `mixed_reduced` reference data already in CI.
- MIXED-species default; add single-species cases only for `cell_caller` and
  count-threshold annotation.
- Enforce content-hash-in-filename from day one and port `check_nf_fixture_hashes.sh`
  (and its `GLOB_SENSITIVE_SUFFIXES`/`REVIEWED_SAFE_SCRIPTS` registries) — external bin scripts
  must be audited for glob/`split('.')` filename parsing before fixtures are hashed.
- Pipeline-stub fixture: a small in-repo `tests/nf/fixtures/*.input_csv.csv` (mirror internal),
  reuse existing `test` profile.

### D. pytest coverage
- Add `pytest.ini` (copy internal: `testpaths = tests`, `--strict-markers --tb=short`,
  markers unit/integration/regression).
- Create `tests/` with per-script test dirs targeting `bin/` logic: `cell_caller.py`,
  `count_matrix.py`, `filter_count_matrix.py`, `categorize_reads.py`, `summary_statistics.py`,
  `features_names.py`, `qc_cascade_plot.py`, `create_single/multi_sample_report.py`,
  `create_curated_barcode_correction_list.py`. Per-dir `conftest.py` to inject `bin/` into
  `sys.path` (external scripts are in top-level `bin/`, simpler than internal's
  `resources/usr/bin/`).
- Use a pixi env under `tests/requirements/` pinning the pipeline's Python/numpy/scanpy/pysam
  (external already ships `pixi.toml` + `pixi.lock` at root — extend or add a test manifest).
- Start with high-value units: cell-calling threshold logic, count-matrix filtering,
  read categorization, summary-statistic arithmetic. Mark heavy scanpy/h5ad tests `integration`.

### E. CI design (extend external `.circleci/config.yml`)
- Add `run-unit-tests` job (pixi -> `pytest tests`), and a `run-nf-tests` job mirroring internal:
  machine executor, `parallelism: N`, `NXF_VER` pinned to 26.04.1, install Java17/Nextflow/nf-test,
  run `check_nf_tests.sh` (warn) + `check_nf_fixture_hashes.sh` (hard fail) on shard 0, then
  `nf-test test --ci --shard ... --junitxml` + `store_test_results`.
- Keep the existing Tower `run-current-branch` job as the e2e tier, but upgrade it to internal's
  pattern: pin `NXF_VER`, poll robustly, and add `tower-validate-s3-outputs` against external's
  real expected outputs (multiqc HTML, summary reports, count matrices) instead of fire-and-forget.
- `check_nf_tests.sh` needs adaptation: internal maps `modules/*/main.nf` -> `modules/*/tests/`;
  external has ONE `modules/processes.nf`, so the mapping must be per-process (grep `^process`
  names) or wait until `processes.nf` is split into per-module dirs during the overhaul.

---

## 6. Portable-to-external (apply without changing external outputs)
- All nf-test scaffolding files, conventions, and the two CI guard scripts.
- Snapshot discipline (text-only snapshots; existence/size for binaries) and stub-test pattern.
- Mixed-species-default fixture strategy + STAR-index-isolation.
- pytest layout (per-module dirs, conftest sys.path injection, markers, pixi-pinned env).
- Tower e2e hardening (NXF_VER pin, status poll, S3 output validation).

## 7. Open questions / decisions for the human
1. Keep external's monolithic `modules/processes.nf`, or split into per-module dirs during the
   overhaul? This determines `testsDir` and how `check_nf_tests.sh` maps changes to tests.
2. External scripts live in top-level `bin/`, not `resources/usr/bin/`. Does the overhaul move
   them under module `resources/` (enabling `moduleBinaries`/staging parity with internal) or
   keep `bin/`? Affects both nf-test config and pytest `sys.path`.
3. Where do external nf-test fixtures live — a new private S3 prefix, or reuse the existing
   public `csgx.public.readonly` data? Fixtures need extraction + content-hashing first.
4. Does external adopt a JS report (and thus Vitest), or keep its current Python HTML report
   (`create_*_sample_report.py`)? Determines whether the JS CI tier is needed.
5. Stub blocks must be ADDED to all 39 processes (currently 0) — is that in scope for the
   overhaul, or a follow-up? Without stubs there is no pipeline DAG test.
6. nf-test snapshots embed an nf-test/Nextflow version in `meta`; regenerate all snapshots under
   the pinned NF 26.04.1 + current nf-test to avoid churn (internal snaps show 25.04.3).
