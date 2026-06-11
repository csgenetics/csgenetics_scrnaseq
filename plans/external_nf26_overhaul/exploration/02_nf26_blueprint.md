# NF26 Migration Blueprint (from INTERNAL pipeline) + External Gap Analysis

Source of truth: `/nssd2/humebc/internal/rnaseq/plans/nf26_migration/` (master_plan.md + phase_1..8, 4b, 4c).
Internal landed: manifest `nextflowVersion = '!>=26.04.0'`, `version = '6.0'` (`internal/rnaseq/nextflow.config:8-9`).
Nextflow installed on beast: 26.04.1 (build 12112).

---

## 1. Categories of change the NF26 migration required

Each with concrete before/after. Phase numbers are the internal plan's.

### A. Quote `env()` arguments in `output:` blocks (Phase 2)
26.04 requires a string-literal argument to `env()`; unquoted was tolerated in 25.x.
- Before: `tuple val(sample_id), env(numreads), path(R1), path(R2)`
- After:  `tuple val(sample_id), env('numreads'), path(R1), path(R2)`
- Internal scope: 18 sites across `gene_annotation/main.nf`, `qc/main.nf`, `cite_seq/main.nf`.
- Verify grep: `grep -rEn "env\([A-Za-z_][A-Za-z0-9_]*\)" modules/` must return zero.

### B. Remove strict-mode / deprecated config flags (Phase 3)
- `nextflow.enable.strict = true` -> DELETE. Its semantics are the default in 26.04. (NB: see Gap Analysis — empirically this flag does NOT error under 26.04.1; internal removed it for tidiness, not because it breaks.)
- `nextflow.enable.configProcessNamesValidation = false` -> DELETE (removed in 26.04). Internal site: `tests/nf/nextflow.config:22`.

### C. Lowercase channel factory (Phase 3)
- Before: `Channel.value([])`  After: `channel.value([])` (canonical lowercase). Applies to `Channel.value/empty/from/of/fromPath/fromList/fromFilePairs`.

### D. Explicit closure params — no implicit `it` (Phase 3)
26.04 requires explicit closure parameters in pipeline operators.
- Before: `.collect { it[1] }`   After: `.collect { it -> it[1] }`
- Internal scope: 16 `.collect { it[N] }` sites in `main.nf`. Convention kept the name `it` (minimal-noise).

### E. `shell:` -> `script:` conversion (Phase 4)
26.04 deprecates the `shell:` directive. Two conversion kinds:
- Group X (full rewrite of `'''...'''` + `!{var}`):
  - Before (shell): `'''  numreads=$(( $(zcat !{R1} | wc -l) / 4))  '''`
  - After (script): `"""  numreads=\$(( \$(zcat ${R1} | wc -l) / 4))  """`
  - Rules: `shell:`->`script:`; `'''`->`"""`; every `!{var}`->`${var}`; every literal bash `$` becomes `\$`; pre-existing `\$` preserved.
  - Result confirmed at `internal/rnaseq/modules/qc/main.nf:37-41` (`count_reads`).
- Group Y (label-rename only): body already `"""..."""` with `${var}`/`\$`; only `shell:`->`script:`.
- Per-process commits (bisect atomicity); regenerate any owned `.snap` inline (only `count_reads.nf.test.snap`).
- Verify: `grep -rEn "^\s*shell:" modules/ main.nf subworkflows/` returns zero.

### F. Process section canonical ordering — `stub:` AFTER executable (Phase 4b)
26.04 strictly enforces section order; `stub:` before `script:`/`shell:`/`exec:` errors with
`Invalid process definition -- check for missing or out-of-order section labels`. Parses fine in 25.x (audit miss).
- Before: `script:` after `stub:`  ->  After: executable block first, `stub:` second. Purely positional; `def` lines and if/else bodies move with the executable block.
- Internal scope: 46 processes / 14 files. Confirmed result `internal/rnaseq/modules/qc/main.nf:13-24` (`merge_lanes`: `script:` then `stub:`).
- Related same-root-cause fix: a `maxForks` directive sitting between `output:` and `script:` errors as `Unrecognized process output qualifier 'maxForks'`; move it to the directive position before `input:`.

### G. `publishDir` directive closure-wrapping (Phase 4c)
26.04 evaluates directive STRING arguments EAGERLY at process compilation (25.x did it lazily at task scheduling). A bare `publishDir "...${sample_id}..."` errors `No such variable: sample_id` because the input var is not in scope at compile time.
- Before: `publishDir "${params.outdir}/.../${sample_id}", mode: 'copy', pattern: '*.tsv'`
- After:  `publishDir { "${params.outdir}/.../${sample_id}" }, mode: 'copy', pattern: '*.tsv'`
- Only the path string is wrapped. A `pattern:` that itself interpolates an input var is also wrapped: `pattern: { "${sample_id}.qc_stats.*" }`.
- `${params.*}`-only strings are SAFE (params resolve at compile time) — do NOT wrap.
- Existing `saveAs: { ... }` are already closures — leave alone.
- Internal scope: 18 directives / 8 files.

### H. Version pin flip + manifest version bump (Phase 5)
- `nextflowVersion = '!>=25.04.0,<26.0.0'` -> `'!>=26.04.0'` (no upper bound; trust semver).
- `version = '5.12'` -> `'6.0'` (major bump signals breaking runtime requirement).
- CI `NXF_VER` bumped to a 26.04.x patch. Single commit (rollback symmetry); the only revert-candidate commit.

### I. CI / Tower launch-layer fixes (Phase 4e — discovered in CI, not syntax)
- Tower auto-selects its platform-default Nextflow (was 25.10.5) which then rejects the `!>=26.04.0` manifest pin. Fix: add `"preRunScript":"export NXF_VER=26.04.1"` to the Tower launch JSON.
- 26.04 checks pipelines out per-revision; `git fetch <sha>` needs the FULL 40-char SHA for non-default-branch commits. Fix: pass `$CIRCLE_SHA1` (full) as the Tower `revision`, keep 7-char SHA only for human-readable name/outdir.

### J. nf-test file fix — closure call sites (Phase 4d)
26.04 v2 strict parser rejects bare `makeDummy(...)` invocation of a `def`-bound closure in workflow-scope test code. Fix: `makeDummy.call(...)`. Test-file only; valid in 25.x too.

### Sequencing principle (portable)
All syntactic changes (A-G) land and stay valid under 25.x BEFORE the version flip (H), so `git bisect` distinguishes migration-syntax regressions from 26.04-runtime regressions. F and G were discovered only at 26.04 PARSE time (`nextflow inspect main.nf -profile test`), not by any 25.x static scan.

---

## 2. What breaks under NF26 strict mode (must fix)

Ordered by how they surface:
1. PARSE-time, surfaced by `nextflow inspect main.nf -profile test` under 26.04:
   - `stub:` before executable block -> `Invalid process definition -- check for missing or out-of-order section labels` (cat F).
   - `maxForks` between `output:` and `script:` -> `Unrecognized process output qualifier 'maxForks'` (cat F).
   - `publishDir`/`pattern:` bare strings interpolating input vars -> `No such variable: <var>` (cat G).
   - Unquoted `env(VAR)` in `output:` (cat A).
   - Capitalised `Channel.value(...)` and implicit-`it` closures (cat C, D).
   - `shell:` directive deprecated (cat E).
   - `nextflow.enable.configProcessNamesValidation` flag removed (cat B).
2. LAUNCH-time (Tower/CI): manifest pin rejection + per-revision full-SHA fetch (cat I).
3. KNOWN-RISK to watch (nextflow#6762): `params is null` for processes reading `params.*` inside script bodies in 26.04 Tower runs. Internal did NOT refactor for it; flagged as halt-and-surface in CI observation.

---

## 3. GAP ANALYSIS — what external STILL needs for full NF26 compliance

External baseline (`/nssd2/humebc/external/csgenetics_scrnaseq`): `nextflowVersion = '!>=25.10.2'`, `version = '1.1.17'` (`nextflow.config:14-15`), `nextflow.enable.strict = true` (`:19`), `process.resourceLimits` present (`:176-180`), no `check_max`. Pipeline is just `main.nf` (473 lines) + `modules/processes.nf` (39 processes).

### Already DONE in external (no action)
- Cat A (env quoting): both `env()` already quoted — `processes.nf:201 env('uniquely_mapped_reads')`, `:519 env('alignment_count')`. Zero unquoted sites.
- Cat C: no capitalised `Channel.` factories; `main.nf:117` already `channel.fromPath`.
- Cat D: all closures already explicit (`.map {row -> ...}`, `{ grouped -> ...}`); zero implicit-`it` operator closures.
- Cat E: zero `shell:` blocks; all 39 processes use `script:`.
- Cat F: zero `stub:` blocks at all -> no section-order or `maxForks` issues.
- `check_max` already removed; `resourceLimits` already configured.

### STILL NEEDED (action required)

GAP-1 (BLOCKER, cat G) — `publishDir` / `pattern:` closure-wrapping.
- EMPIRICALLY CONFIRMED: `NXF_VER=26.04.1 nextflow inspect main.nf -profile test` FAILS:
  `ERROR ~ No such variable: sample_id -- ... modules/processes.nf at line: 498`.
- Sites needing the path string wrapped in `{ ... }` (interpolate `${sample_id}`): 14 publishDir path strings in `modules/processes.nf` at lines 533, 534, 535, 677, 678, 679, 736, 737, 738, 739, 740, 793, 823, 866.
- Sites needing the `pattern:` argument wrapped (interpolate `${sample_id}`): 2 — line 498 (`pattern: "${sample_id}.annotated.bam"`) and line 711 (`pattern: "${sample_id}*_pdf_with_cutoff.html"`).
- LEAVE ALONE: 10 `${params.outdir}`-only path strings (compile-time safe). Existing `saveAs: {...}` closures at lines 269, 535, 562 are already closures.
- NOTE line 269: path is params-only but `saveAs` interpolates `${sample_id}.${prefix}` inside an already-closure -> safe, no change.

GAP-2 (cat B) — remove `nextflow.enable.strict = true` (`nextflow.config:19`).
- Strict is the default in 26.04, so the flag is redundant. EMPIRICAL CAVEAT: under 26.04.1 the flag does NOT currently error (`nextflow config -profile test` exits 0 with `strict = true` echoed). So this is cleanup/idiomatic, not a hard blocker. Internal deleted it; external should follow for parity. External has NO `configProcessNamesValidation` to remove.

GAP-3 (cat H) — version pin + manifest bump.
- `nextflow.config:14`: `'!>=25.10.2'` -> `'!>=26.04.0'`.
- `nextflow.config:15`: `version` bump. Decide major-vs-minor for external's customer-facing `1.1.17` (internal chose a MAJOR bump to signal the breaking runtime requirement; external's semantics are customer-facing so the bump policy is a human decision — see open questions).
- Update any CI `NXF_VER` pin (external CI config not yet located in this pass — must be found and bumped; internal used `.circleci/config.yml:148`).

GAP-4 (cat I, CI/Tower) — if external launches via Tower/Seqera, apply the same preRunScript `export NXF_VER=26.04.x` and full-40-char-SHA revision fixes. Needs verification of external's CI/launch mechanism.

GAP-5 (cat J / nf-test) — external's nf-test suite (if any) must be re-run under 26.04 and snapshots (`meta.nextflow` field) regenerated; watch for the closure-call-site pattern. External test layout not assessed in this pass.

### After GAP-1..3, re-verify
`NXF_VER=26.04.1 nextflow inspect main.nf -profile test` must parse clean. Then full test sweep on 26.04, watching for nextflow#6762 (`params is null`).
