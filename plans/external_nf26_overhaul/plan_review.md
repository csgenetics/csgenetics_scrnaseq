# Consolidated plan review — external NF26 overhaul

Lead-reviewer consolidation of five adversarial reviews. Each finding below was vetted against the
actual repository files (`modules/processes.nf`, `main.nf`, `bin/assign_multi_mappers.gawk`,
`bin/create_multi_sample_report.py`, `.circleci/config.yml`, `nextflow.config`) and the three plan docs.
Overlapping findings are merged; unsubstantiated ones are listed at the end with the reason.

Severity counts: **5 blockers, 6 majors, 4 minors.**

---

## Blockers

### B1. Phase 0 captures the golden baseline on the wrong runtime — it cannot run as written
Phase 0 (`master_plan.md:41`) says "pin NF26, run `-profile test` on CURRENT code". But the plan's own
confirmed fact (`master_plan.md:19-20`) is that current code FAILS to parse under NF26
(`No such variable: sample_id` at `processes.nf:498`), and `regression_strategy.md:10-11` correctly
states the baseline must be produced on Nextflow **25.10.2**. As written, Phase 0 halts at the first
`inspect`/run and no baseline is ever captured, blocking every downstream output-equivalence check.
**Fix:** Rewrite the Phase 0 row to defer to `regression_strategy.md §1` verbatim: tag `baseline-pre-nf26`;
run `test` + `test_hsap_mmus_2_sample_full` + the degenerate fixture on CURRENT code under **25.10.2**;
run twice for the determinism probe (§5); archive `results/` + `manifest.sha256` + `trace.txt`.

### B2. STAR `--runThreadN` is mis-classified R0 — multimapper representative selection is order-sensitive
`optimization_assessment.md:27` classifies `--runThreadN 8 -> ${task.cpus}` as R0 ("STAR BAM is
thread-count-independent"). Verified: `processes.nf:204-213` runs STAR with NO
`--outMultimapperOrder`/`--outSAMmultNmax` override, so multimapping-alignment order is STAR-default and
thread-dependent. That order propagates through `processes.nf:410/451` (`samtools sort -n | gawk`) into
`assign_multi_mappers.gawk`. Verified in the gawk: `assigned_array_corrected[$18]` is keyed ONLY by the
XT gene tag (`bin/assign_multi_mappers.gawk:72`), and the `END` block emits to `assigned_reads.sam_body`
only when `length(assigned_array_corrected)==1` (line 94). When a read has two "Assigned" alignments to
the SAME gene, last-write-wins both collapses the map to length 1 (flipping the read from ambiguous to
assigned) AND picks whichever alignment arrived last — so the surviving record's coordinates depend on
input order, which depends on thread count.
**Fix:** Reclassify STAR threading as **R1** (must-be-proven). Either pin determinism explicitly
(`--outMultimapperOrder Random --runRNGseed <fixed>`) and prove the gawk selection is order-invariant,
or make the determinism probe vary thread count (baseline at 8 vs `${task.cpus}`) and compare the
DOWNSTREAM count-bearing outputs (`bcGeneSummary.txt` multiset, h5ad `X`), not just canonically-sorted
BAMs. Add an nf-test feeding the gawk the same read group in two input orders.

### B3. umi_tools dedup representative selection is order-sensitive — the threads-fixed probe cannot catch it
Verified: `processes.nf:625-628` runs `umi_tools dedup --per-cell` picking one representative alignment per
(cell, position, UMI-network); its XT tag is extracted by `io_count` (`processes.nf:661`) into the
published `bcGeneSummary.txt` -> h5ad `X` (Class C, compared EXACT per `regression_strategy.md:31`). The
directional network is deterministic on the UMI multiset, but the emitted representative is selected by
input traversal order. The determinism probe (`regression_strategy.md:58-63`) re-runs with IDENTICAL
threads, so input order to dedup is byte-identical between probe runs; it never exercises the
order-perturbed path that the Tier-1 threading/`samtools sort` changes introduce. Class-B BAM
normalization absorbs the BAM record-order difference but NOT its count consequence.
**Fix:** Add a dedicated order-perturbation test independent of the fixed-thread probe: feed one fixture
BAM to dedup in two record orders and assert `bcGeneSummary.txt` and h5ad `X` are identical. If not, the
Tier-1 threading changes are NOT R0 and must be gated behind a deterministic-ordering fix.

### B4. moduleBinaries split breaks the two `file("${baseDir}/bin/assign_multi_mappers.gawk")` references
Verified: `main.nf:304` and `main.nf:306` pass the gawk as an explicit `file()` process INPUT, received as
`path(multi_mapper_script)` and invoked `gawk -f $multi_mapper_script` (`processes.nf:410,451`).
`master_plan.md:28,43` move `bin/` under module `resources/usr/bin/` with `nextflow.enable.moduleBinaries`.
moduleBinaries only injects the owning module's `resources/usr/bin` onto PATH; it does NOT make a file
resolvable via `file("${baseDir}/bin/...")`. The plan's GAP analysis never inventories
`${baseDir}`-relative asset references (also templates at `main.nf:179,180,182`).
**Fix:** Keep `assign_multi_mappers.gawk` in top-level `bin/` (it is shared by two processes and consumed by
`file()`, not PATH). If it must move, rewrite both `file()` paths in `main.nf` to the new location and
re-run `-profile test` to confirm byte-identical multimapper BAMs. Before any move, grep `main.nf` for all
`${baseDir}`/`${projectDir}` references and enumerate them in the blueprint.

### B5. moduleBinaries split breaks cross-script Python imports in `create_multi_sample_report.py`
Verified: `bin/create_multi_sample_report.py:23` does
`from create_single_sample_report import get_cell_stat_cat_dict_obj, read_html_plot` (and imports from
`cell_caller`). The `multi_sample_report` process (`processes.nf:885,901`) invokes the script bare-name,
relying on Nextflow's implicit whole-`bin/`-on-PATH staging. But `create_single_sample_report.py`
(single_summary_report, `processes.nf:863`) and `cell_caller.py` (cell_caller, `processes.nf:705`) belong
to DIFFERENT processes and would land in different module `resources/usr/bin` dirs under moduleBinaries —
so the imports resolve to nothing and the report process fails at runtime (ImportError). Nothing tests
`multi_sample_report` before Phase 4.
**Fix:** Inventory all intra-`bin/` imports first. Either keep these scripts co-resident in one shared
module's `resources/usr/bin`, package shared logic into a real importable module staged into every
consumer, or do NOT move `bin/` under moduleBinaries (keep top-level `bin/`). Add an nf-test for
`multi_sample_report` in the SAME wave as the split, not Wave 4.

---

## Majors

### M1. Phases 2 and 3 are declared parallel but both mutate `modules/processes.nf` — guaranteed collision
`master_plan.md:55` claims Phases 2/3/5 are independent. But Phase 2 (`master_plan.md:43`) splits/deletes
`processes.nf`, while Phase 3 (`master_plan.md:44`) edits line-level `script:`/`cpus` inside the same
processes (star `processes.nf:205`, `sort_index_bam` ~599, `initial_feature_count` 307/310) AND
`conf/base.config`. A whole-file move concurrent with line edits to the same file is unmergeable.
**Fix:** Serialize: Phase 2 (structural split, output no-op) lands FIRST, then Phase 3 edits the
already-split per-process module files. Drop the claim that 2 and 3 are independent.

### M2. Report consolidation (Decision 2 Option C) contradicts the comparator's published-file-SET contract
`master_plan.md:30` / `exploration/05_report.md:104` Option C folds/removes published HTML
(`multisample_summary_plots.html`, `multisample_qc_cascade.html`, per-sample `${sample_id}_report.html` at
`processes.nf:873,896-897`). `regression_strategy.md:51-52` FAILS loudly on any missing/extra published
path. The "audit downstream consumers" escape hatch is hollow for an EXTERNAL pipeline — consumers are
customers' private scripts/LIMS that the team cannot enumerate, so the audit never returns "clear", yet
the decision is LOCKED.
**Fix:** Resolve before Phase 5: either (a) constrain Option C to in-place restyle that keeps the EXACT
published filename set (forbid folding/removing files), or (b) treat removed HTML as a deliberately-blessed
output change with its own re-baseline commit + customer deprecation note (same gating as the dedup swap).
Do not leave "preserve OR audit" as a free mid-flight choice.

### M3. The gated Tier-3 / dedup re-baseline has no phase and collides with Phase 5's metric freeze
`master_plan.md:32-34` puts the Rust fusions + `umi_tools dedup -> count_table_builder` swap in scope, but
the only Phase 3 row (`master_plan.md:44`) lists Tier 1 R0 changes only; `optimization_assessment.md:67-70`
calls Tier 3 a separate gated workstream with no phase number. The dedup swap changes
h5ad/called-cells/dedup.log -> `metrics.csv`. Phase 5 (`master_plan.md:46`) runs in parallel and freezes
golden metrics; if Phase 5 freezes while a blessed dedup re-baseline lands, the golden files are stale at cut.
**Fix:** Add an explicit late phase (e.g. 3b) for Tier 3, depending on Phase 4 nf-test gates and user
sign-off (`regression_strategy.md §8`), sequenced AFTER Phase 5's freeze, re-cutting only the
dedup-downstream golden files. State that Phase 3 ships Tier 1 only.

### M4. bcGeneSummary.txt is Class A "exact line compare" but its line order is volatile dedup.bam order
`regression_strategy.md:27` classes `*_bcGeneSummary.txt` as Class A exact, and §3 normalization
(lines 38-39) does NOT sort lines. Verified: `processes.nf:661` emits one line per XT-tagged dedup record
in dedup.bam ORDER with no post-sort, and the plan itself admits that order is volatile (line 63). An exact
line compare FAILS on benign reordering (phantom regression); naive line-sorting would then mask a real
representative-swap that changes a gene while keeping the same count.
**Fix:** Reclassify with explicit normalization: compare the SORTED multiset of (io_sequence, gene_name)
tuples (the count map), not the raw line stream — order ignored, gene swaps still caught.

### M5. nf-test + fixture-hash CI needs AWS creds / `machine` + `large` resource_class that external CircleCI lacks
Verified: external `.circleci/config.yml` has ONE `docker` job `run-current-branch` on `resource_class:
small` whose only secrets are `TOWER_*`. Internal's nf-test job is `machine`, `resource_class: large`,
`parallelism: 4` with S3 fixtures needing `AWS_ACCESS_KEY_ID`. External `conf/aws.config` hardwires
`awsbatch` which the nf-test config must override. `master_plan.md:45` Phase 4 / Phase 6 never list these.
**Fix:** Add a Phase 6 prerequisite: provision AWS read creds as CircleCI project env vars (or confirm
fixtures live in a public read-only bucket with unsigned S3 so no creds are needed), add
`machine`+`resource_class: large`+`parallelism` to the new job, and ensure `tests/nf/nextflow.config` does
NOT `includeConfig conf/aws.config`.

### M6. Per-PR CI gates only single-species `test`; the mixed-species path has no per-PR-runnable golden
`regression_strategy.md §9` runs the comparator per-PR for the `test` profile only; mixed runs "on a
schedule or label". The only mixed coverage is `test_hsap_mmus_2_sample_full`, which points at the full
production-scale mixed STAR index — too heavy for per-PR. Mixed-only logic (`cell_caller --single_species
false` `processes.nf:724`, `mixed_args` in count/filter, species-tagged GTF) is exercised by no PR check,
so a module-split/threading drift in the barnyard path passes CI and ships.
**Fix:** Build a REDUCED mixed-species test profile (chr-subset human+mouse index, 2 tiny samples) sized to
run per-PR, capture its golden, and gate it alongside the single-species golden.

---

## Minors

### m1. No phase covers customer docs/README, CHANGELOG, version-bump, or rollback for a BREAKING NF26 floor
GAP-3 flips `nextflowVersion '!>=25.10.2'` (verified `nextflow.config:14`) to `!>=26.04.0`, forcing all
customers onto NF26, and the manifest `version = '1.1.17'` (verified `nextflow.config:15`) bump is still an
open human decision. `master_plan.md` has no readme/changelog/rollback/version-policy phase; Phase 10 only
opens a PR.
**Fix:** Add a docs/release phase before Phase 10: update README runtime requirement, add a CHANGELOG entry
stating the breaking floor and resolving the 1.1.17 major-vs-minor question, and document the rollback
(the GAP-3 version-pin flip is the revert-candidate commit).

### m2. dedup.log Class A normalization mishandles umi_tools' non-ISO timestamp and the empty-branch literal `\n`
`regression_strategy.md:27` classes `*.dedup.log` Class A. umi_tools emits comma-millisecond timestamps
(not ISO-8601) with counts on the same INFO line, so a naive line-strip risks deleting the count line.
Verified: the empty branch `processes.nf:633` writes
`echo "INFO Reads: Input Reads: 0\nINFO Number of reads out: 0\n"` WITHOUT `-e`, so `\n` stays literal —
one physical line, a different shape from the real umi_tools log.
**Fix:** Normalize dedup.log to ONLY the two integers `summary_statistics.py` consumes (input reads, reads
out); discard other lines. Verify the empty-branch literal-`\n` output parses to the same two integers.

### m3. check_nf_tests.sh modified-process->test mapping is undefined on the monolithic processes.nf
Internal maps `modules/*/main.nf` -> `modules/*/tests/`. External has ONE `processes.nf` with ~33-39
processes. The guard cannot map "modified process needs a test" on the monolith, but Phase 3 can edit
processes before Phase 4's guard exists.
**Fix:** Land the Phase 2 split strictly before any process-editing phase, OR write an interim guard that
greps `^process <name>` from the `processes.nf` diff and requires a matching test file.

### m4. GAP-1 publishDir site list labelled "EMPIRICALLY CONFIRMED" but inspect only validated line 498
`nextflow inspect` halts at the first error (498), so the other 15 sites are static-audit, not
inspect-validated. (The list is correct, but the provenance claim is overstated.)
**Fix:** Downgrade wording to "one site confirmed by inspect; remainder by static audit"; in Phase 1, wrap
and re-run inspect iteratively until clean.

---

## Discarded / downgraded findings

- **Seqera Phase 7 CE/credentials/role/bucket findings (4 findings: external-CE creation, pipelineuser
  launch-only role, TOWER_COMPUTE_ENV_ID rewiring, branch-name vs full-SHA revision).** Not discarded as
  wrong — the CircleCI evidence is verified (branch-name `revision`, `awk -F:` parsing,
  `s3://csg-tower-bucket`/`csg-nextflow` hardcoded, single small job). They are genuine. They are folded
  into the plan as a single consolidated caveat rather than four separate issues, because they all reduce
  to one actionable amendment: **Phase 7 must not assume "mirror internal" is copy-paste.** Treated as
  major-equivalent operational risk, but they are external-credential/infra prerequisites the agent cannot
  self-serve and do not block the code overhaul; captured under M-class amendments to the Phase 6/7 rows
  (see amendments). Listed here to avoid double-counting in the severity totals while still requiring the
  edit.
- **quay.io push-access finding.** Verified-plausible (login != push scope) but unprovable from the repo;
  folded into the same Phase 7/8 credential-verification amendment rather than counted separately.
- **resolved_configuration.txt cross-runtime JSON formatting drift (note).** Plausible but the reviewer
  could not diff 25.10.2 vs 26.04.1 (only 26.04.1 on beast). Real risk; folded into m-class as a
  determinism-probe amendment rather than a standalone issue.
- **Class D Plotly numeric re-extraction is redundant (note).** Substantiated as a design observation
  (Jinja report values derive from `metrics.csv` already in Class A; KDE traces derive from h5ad already
  in Class C). Not a defect that ships bad output — it is a brittleness/effort note. Downgraded to an
  amendment on `regression_strategy.md §3` Class D, not counted as an issue.
- **Degenerate fixture conflates two empty-trigger branches.** Verified: `main.nf:243` (`countFastq()==0`,
  post-QC) and `main.nf:256` (STAR `uniquely_mapped_reads==0`) are two independent feeds and
  `processes.nf:305` uses lexical `[[ $aligned_count > 0 ]]`. Real, but it is a fixture-construction detail
  inside §6, not plan-blocking; folded into the M1/B1 baseline amendment (split into two fixtures).

---

## Plan amendments required

### master_plan.md
1. **Phase 0 row (line 41):** rewrite to defer to `regression_strategy.md §1` verbatim — baseline on
   **Nextflow 25.10.2** (not NF26), run `test` + `test_hsap_mmus_2_sample_full` + degenerate fixture, run
   twice for the determinism probe, archive `results/` + `manifest.sha256` + `trace.txt`. (B1)
2. **Parallelism (line 55):** drop the "Phases 2,3 independent" claim; serialize Phase 2 (split) before
   Phase 3 (threading). (M1)
3. **Decision 1 / Phase 2 (lines 28,43):** do NOT move `bin/assign_multi_mappers.gawk` (file()-consumed by
   2 processes) or the cross-importing report scripts (`create_multi_sample_report.py` ->
   `create_single_sample_report`/`cell_caller`) under moduleBinaries; keep them top-level, or rewrite all
   `file()` refs and co-resident imports and add an nf-test in the split wave. (B4, B5)
4. **Decision 2 (line 30):** replace "preserve filenames OR audit consumers" with a hard choice — in-place
   restyle keeping the exact published-file SET, or a blessed re-baseline with deprecation note. (M2)
5. **Add Phase 3b:** gated Tier-3 / dedup re-baseline, dependent on Phase 4 + user sign-off, sequenced
   AFTER the Phase 5 metric freeze; Phase 3 ships Tier 1 only. (M3)
6. **Phase 6 row + Credential column:** add AWS read creds + `machine`/`resource_class: large`/parallelism
   for nf-tests, and exclude `conf/aws.config` from the nf-test config chain. (M5)
7. **Phase 6/7/8 rows:** add explicit pre-flight checks — pipelineuser role (CE-create/credentials-register
   vs launch-only), register an external-workspace AWS-credentials object + work bucket before CE creation,
   out-of-band CircleCI env-var update (`TOWER_COMPUTE_ENV_ID`, confirm `TOWER_WORKSPACE_ID`), quay push
   scope probe, and copy internal's full-SHA `revision` + `NXF_VER` preRunScript fixes. Mark these as
   human/external prerequisites, not "mirror internal" copy-paste. (Discarded-but-folded Seqera/quay items)
8. **Add a docs/release phase before Phase 10:** README runtime bump, CHANGELOG breaking-floor entry +
   version decision, documented rollback. (m1)

### optimization_assessment.md
9. **Tier 1 table, `star` row (line 27):** reclassify `--runThreadN` from R0 to **R1**; remove the
   "STAR BAM is thread-count-independent" assertion; require either pinned multimapper ordering or a
   thread-varying determinism probe comparing downstream counts. (B2)
10. **Tier 1 verification note (lines 36-37):** state that the determinism probe with FIXED threads cannot
    prove R0 for any change that perturbs BAM record order feeding the gawk or umi_tools dedup; add the
    order-perturbation test as the gate. (B2, B3)

### regression_strategy.md
11. **§2/§3 Class A:** reclassify `bcGeneSummary.txt` to compare the sorted (io_sequence, gene_name)
    multiset, not the raw line stream. (M4)
12. **§3 Class A dedup.log:** normalize to the two integers `summary_statistics.py` consumes; verify the
    empty-branch literal-`\n` output parses to the same two integers. (m2)
13. **§5 determinism probe:** add an order-perturbation test independent of the fixed-thread probe (feed
    dedup/gawk the same records in two orders; assert identical `bcGeneSummary.txt` + h5ad `X`). Also render
    `resolved_configuration.txt` under both 25.10.2 and 26.04.1 and add a JSON-key-map normalization if the
    serialization differs. (B2, B3, resolved_config note)
14. **§6 branch-coverage matrix:** split the degenerate row into two fixtures — a QC-empties-to-zero-reads
    FASTQ (`main.nf:243`) and a survives-QC-but-zero-unique-alignments FASTQ (`main.nf:256`); note
    `processes.nf:305` uses lexical `>` so the fixture must drive `aligned_count` to literal `0`. Add a
    reduced per-PR mixed-species profile + golden. (B1, M6)
15. **§3 Class D:** for Jinja reports assert structural presence (sample set, metric KEYS) only — value
    re-extraction is redundant with Class A `metrics.csv`; for Plotly compare only non-derivable scalar
    annotations, not full KDE trace arrays. (Class-D note)
