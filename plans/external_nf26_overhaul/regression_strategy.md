# Output-equivalence regression strategy

> AMENDED by adversarial review — apply `plan_review.md` amendments 11-15 before using this doc as the build
> spec: (11) `bcGeneSummary.txt` compared as a SORTED (io_sequence, gene_name) multiset, not raw line order;
> (12) `dedup.log` normalized to the two integers `summary_statistics.py` consumes; (13) §5 probe must ALSO
> perturb record order (not just rerun with fixed threads) to expose order-sensitive count drift, and render
> `resolved_configuration.txt` under both 25.10.2 and 26.04.1; (14) §6 split the degenerate row into TWO
> fixtures (QC-zero-reads `main.nf:243` vs survives-QC-zero-unique-alignments `main.nf:256`) and add a REDUCED
> per-PR mixed-species profile + golden; (15) §3 Class D asserts structural presence only (values are already
> covered by Class A `metrics.csv` / Class C h5ad).


Purpose: guarantee that every rewrite (NF26 migration, module split, Tier 1-2 optimizations, and any
adopted Tier 3 fusion) produces outputs that are identical to a frozen baseline — or, where a change is
deliberately output-altering (the dedup engine swap), that the change is quantified and explicitly blessed.

This is the linchpin of the whole overhaul. Nothing in Phase 3+ lands without passing it.

## §1. Frozen baseline (golden reference)
- **Code point:** tag the pre-change commit (current `devel` head) as `baseline-pre-nf26`. The baseline is
  produced from THIS code, run on Nextflow **25.10.2** (the version the current code parses under).
- **Profiles run for the baseline** (coverage rationale in §6):
  - `test` — fast human chr1-subset, 2 samples (smoke + the default golden set).
  - `test_hsap_mmus_2_sample_full` — the ONLY mixed-species/barnyard path.
  - a small **degenerate-sample** fixture (≥1 sample with 0 reads / 0 alignments) to cover the
    `create_valid_empty_bam` + `aligned_count==0` branches (build if not already present).
- **Archive per profile:** the complete published `results/` tree, plus:
  - `manifest.sha256` — sha256 of every published file (raw).
  - `manifest.normalized.json` — per-file normalized fingerprint (see §3) for files with volatile bytes.
  - `trace.txt` (Nextflow trace) for per-process wall-time/CPU/RAM — this also ranks the real hotspots.
- Store baselines out-of-tree (e.g. `/nssd2/humebc/external-nf26-runs/baseline/<profile>/`) and commit
  ONLY the small text golden files + manifests under `tests/regression/golden/<profile>/`.

## §2. Output classification (drives the comparator)
Every published output is assigned a class once, in a registry `tests/regression/output_classes.yml`:
- **A — deterministic text:** `*.metrics.csv`, `multisample_out.csv`, `*_bcGeneSummary.txt`,
  `resolved_configuration.txt`, RSeQC `*_RSeQC.txt`, `*.dedup.log`, `barcodes/features.tsv.gz`,
  `matrix.mtx.gz`. → exact compare after §3 normalization. Counts/integers EXACT.
- **B — deterministic content, volatile bytes:** all BAMs. Record SET and counts are deterministic given
  identical input; byte order / `@PG` / compression are not. → samtools-normalized compare (§3).
- **C — binary structured:** `*.h5ad` (anndata). → compare `X` (sparse counts) EXACT, `obs`/`var` indices
  EXACT, `layers`/`uns` value-equal; ignore HDF5 container byte layout.
- **D — plots / reports:** `*.html` (cell-caller, qc_cascade, multiqc, summary reports), `*.png`. → NEVER
  byte-compare. Extract the underlying numeric series / embedded metric values and compare those; assert
  structure (same plot divs, same sample set) exists.

## §3. Normalization rules (so we compare signal, not noise)
- **Class A:** strip/booleanize volatile fields — absolute work/output paths, ISO timestamps, tool version
  strings, Nextflow session UUIDs, hostnames. `resolved_configuration.txt` normalizes `outdir`/`workDir`.
- **Class B (BAM):** `samtools sort` to a canonical order, `samtools view` stripping `@PG`/`@CO`/`@HD SO`,
  then sha256 of the record stream. Plus independent invariants that must match EXACTLY:
  `samtools flagstat`, `samtools idxstats`, `samtools view -c`, per-tag counts used downstream
  (`[NH]`, `[nM]`, `[XS]`, `XT` presence). For Tier-1 threading changes these are byte-stable modulo `@PG`.
- **Class C (h5ad):** load with anndata, compare `adata.X` via sparse equality (same nnz, same indices,
  same data), `var_names`/`obs_names` set+order, and any metrics-bearing `obs`/`var` columns.
- **Class D:** for Plotly HTML, parse the embedded JSON `data` traces and compare numeric arrays; for the
  Jinja reports, parse out the rendered metric values and assert they equal the matching `metrics.csv`.

## §4. The comparator harness
`tests/regression/compare_outputs.py <baseline_dir> <candidate_dir>`:
1. Pairs files by relative path under `results/`; FAILS loudly on any missing/extra published path
   (the published-file SET is itself part of the contract).
2. Dispatches each pair to its class comparator (§2/§3).
3. Emits `regression_report.json` + a human summary: per-file PASS/FAIL with the specific diff
   (which metric, which gene, which count) on failure.
4. Exit non-zero on any FAIL. No silent tolerance — every non-exact comparison is an explicit, named rule.

## §5. Determinism probe (establish the noise floor FIRST)
Before trusting any "EXACT" claim: run the **baseline twice** (same code, same NF 25.10.2, same inputs) and
run the comparator baseline-vs-baseline. Anything that differs is inherently non-deterministic and its
normalization rule (§3) must absorb that difference — otherwise we'd chase phantom regressions. This
calibrates the comparator and is a prerequisite to Phase 3. Expected non-determinism: `@PG` command lines,
timestamps, plot div ids, possibly multi-threaded BAM record order (which is exactly why class B normalizes).

## §6. Branch-coverage matrix (the baseline must exercise what rewrites touch)
| Code branch | Covered by |
|---|---|
| single-species (human) full path | `test`, `test_pbmc_4_sample_full` |
| mixed-species / barnyard (`mixed_species=true`) | `test_hsap_mmus_2_sample_full` |
| empty / 0-alignment sample (`create_valid_empty_bam`, `aligned_count==0` branches) | degenerate fixture (build) |
| multimapper present + absent (`assigned_reads.sam_body` exists / not) | covered by real data in `test` + mixed |
| manual cell-caller threshold path | add a small fixture using `manual_threshold_template.csv` |
A rewrite that touches a branch with no baseline coverage MUST get a fixture before it lands.

## §7. Two-level verification per change
- **End-to-end:** re-run the affected profile(s) on the candidate code (NF 26.04.1 from Phase 1 onward),
  run §4 comparator vs the frozen baseline. Must PASS.
- **Per-process (nf-test):** every rewritten process gets an nf-test that runs it on a committed fixture and
  snapshots the normalized output (class-appropriate: text snapshot for A, flagstat/idxstats/record-hash for
  B, X-matrix summary for C). For Tier 1 this localizes any drift; for Tier 3 it is the equivalence gate
  comparing the new binary's output against the OLD chain's output on the same fixture.

## §8. Special case — deliberately output-changing dedup (Tier 3 `count_table_builder`)
The `umi_tools` → Rust UMI-grouper swap is expected to change outputs. It is NOT verified by identity but by
a quantified, bounded, blessed delta:
- Per sample, report: Δ total deduplicated reads, Δ per-cell UMI distribution (and its correlation),
  Δ genes-detected-per-cell, Δ number of called cells, and a per-gene count scatter/correlation.
- Define acceptance bounds WITH the user before adoption (e.g. called-cells within ±X%, per-gene count
  Pearson r > 0.999). Present the numbers; only adopt on explicit sign-off, and re-cut the golden baseline
  for the dedup-downstream outputs at that point (recorded as an intentional re-baseline commit).

## §9. CI gating
- A `regression` CI job runs the comparator against the committed golden manifests for the `test` profile on
  every PR (full/mixed profiles run on a schedule or label, being heavier).
- Golden files are regenerated only via an explicit `--update-baseline` path that is reviewed in the PR diff;
  CI never auto-updates them. A changed golden file in a PR is a red flag requiring justification.

## §10. What "done" means for a rewrite
A rewrite is accepted only when: (a) determinism probe calibrated; (b) end-to-end comparator PASS on every
covering profile; (c) its per-process nf-test PASS; (d) for R2, the equivalence verdict is EXACT, or the
delta is quantified and user-blessed. Anything else halts and is surfaced.
