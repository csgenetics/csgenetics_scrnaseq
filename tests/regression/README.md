# Output-equivalence regression comparator

`compare_outputs.py` decides whether the published-output trees of two pipeline
runs are *output-equivalent*. It is foundational verification infrastructure for
a customer-facing pipeline: it fails loud and applies **no silent tolerance**.
A real divergence in any class that is meant to be exact drives a nonzero exit
code.

## Usage

```bash
python3 compare_outputs.py <dir_a> <dir_b> [--json OUT] [--max-diffs N] \
                           [--envelope-max-flips N] [--allow-subset]
```

- `<dir_a>`, `<dir_b>` — the two output directories to compare.
- `--json OUT` — also write the full per-file result set (with the run summary)
  as JSON to `OUT`.
- `--max-diffs N` — non-negative cap on the number of differing
  lines/triplets/keys reported per file (default 10; zero suppresses those
  details). It only bounds the *report*, never the verdict.
- `--envelope-max-flips N` — enable **envelope mode** for the count matrices
  (`MTX`, `H5AD`). Unset (default) = strict mode. See
  [Envelope mode](#envelope-mode-the-count-matrix-ambiguity).
- `--allow-subset` — explicitly compare a curated, non-empty common subset.
  This supports diagnostic 1.x-to-2.0 inventories, where intentional filename
  changes make whole-tree path equality impossible; it does not relax any
  selected file's class. The default compares two 2.0 published trees and
  requires all five fixed `pipeline_info` files.

Exit code is `0` when equivalent, `1` when any failure is detected.

Dependencies: Python standard library only for output trees without H5AD or BAM
files. Comparing `.h5ad` requires `anndata`; comparing `.bam` requires
`samtools` on `PATH`. If a required reader is unavailable, the affected file is
reported `DIFFER` and the comparison fails. Skipping validation would allow
identical corrupt bytes to pass.

CI runs the adversarial comparator suite with the committed Pixi lock at
`tests/requirements/comparator/pixi.lock`, which pins Python, pytest, anndata,
numpy, pandas, scipy, and samtools. To reproduce that gate locally:

```bash
pixi run --frozen --manifest-path tests/requirements/comparator/pixi.toml \
  python -m pytest tests/python/test_compare_outputs.py -v
```

## Contract: the published-file set

The *set* of relative paths is part of the contract. If any relative path exists
under one directory but not the other, that is a `FILE_SET_MISMATCH` and the run
fails. Extra or missing published files are never tolerated.

By default, `resolved_configuration.txt`, execution trace, report, timeline and
DAG must all exist beneath `pipeline_info/` in both trees. An empty tree, or two
trees both missing any of those fixed manifest files, therefore fails.
`--allow-subset` disables only that fixed-manifest requirement; the curated
subset must still be non-empty and its relative file set must match exactly.
Dynamic sample outputs cannot be inferred without the run inputs, so they are
governed by symmetric path-set equality rather than a hard-coded minimum list.

## Per-class rules

Each paired file is classified by the **first** matching rule below, then
compared by that class's logic. Only explicitly documented run-specific,
benign volatility is normalised away. “Strict mode” means that every file must
satisfy this per-class contract; it does **not** mean that every file in the
output tree must have identical container bytes.

| Class | Matches | Comparison | DIFFER = failure? |
| --- | --- | --- | --- |
| `TEXT_EXACT` | `*.metrics.csv`, `multisample_out.csv`, `*.txt` under an `RSeQC/` path | require UTF-8, then exact byte compare; reports up to `--max-diffs` differing lines | yes |
| `DEDUP_LOG` | `*.dedup.log` | parse either the exact current two-line read-count summary or the explicit completed UMI-tools 1.1.2 legacy shape; validate `0 <= reads out <= input reads` (and the legacy positions count), then compare input/output tuples | yes |
| `CONFIG` | `resolved_configuration.txt` | parse a non-empty parameter object, require valid `outdir` and stable payload, validate optional `workDir`, then drop those run paths and compare the remainder | yes |
| `TRACE` | `execution_trace.txt` | parse the pinned Nextflow 26.04.1 TSV schema; compare the multiset of process/task/tag/status/exit records; reject malformed, incomplete, failed, or nonzero-exit traces | yes |
| `MULTIQC_LOG` | `multiqc.log` | parse the pinned MultiQC 1.14 log; normalise timestamps, run/temp directories and the external update check; compare all remaining ordered log records and reject error-level/incomplete logs | yes |
| `MULTIQC_JSON` | `multiqc_data.json`, `*.multiqc.data.json` | strictly parse the MultiQC 1.14 data object; normalise creation/analysis paths and source directories; compare all report data/configuration that remains | yes |
| `MULTIQC_SOURCES` | `multiqc_sources.txt` | parse the four-column TSV and compare the multiset of module/section/sample/source-basename records | yes |
| `GZ_TEXT_EXACT` | `*barcodes.tsv.gz`, `*features.tsv.gz` | gunzip, exact compare | yes |
| `MTX` | `*matrix.mtx.gz` | strictly validate MatrixMarket structure and exact non-negative numeric values, then compare shape/field and entry triplets as an externally sorted multiset | yes |
| `H5AD` | `*.h5ad` | open with anndata; require finite, non-negative real-numeric `X` and `raw.X`; compare `X`, obs/var annotations and names, layers, multidimensional annotations, pairwise matrices, `uns`, and `raw`; named `*.empty.h5ad` sentinels pass only when both are zero bytes | yes |
| `BAM` | `*.bam` | fully read with samtools; compare the stable SAM header and the multiset of molecule identities within each coordinate tie | yes |
| `HTML` | `*.html` | parse elements/scripts and apply a filename-specific consolidated-report, MultiQC, Plotly, or Nextflow contract; compare stable Nextflow task/timeline/DAG semantics; explicitly named empty Cell Caller plots pass only when both are zero bytes | invalid, failed, unrecognised, or semantically changed Nextflow HTML fails; valid presentation HTML is `PRESENT` |
| `BINARY_EXACT` | anything else | sha256 compare | yes |

### Why these normalisations

- **`DEDUP_LOG`** — the contig-parallel producer consolidates volatile umi_tools
  logs into the two fields consumed by summary statistics: `Input Reads` and
  `Number of reads out`. Both must occur exactly once, and the output count
  cannot exceed the input count. The empty-sample branch emits the same two-line
  schema with zero counts. For the documented 1.x-to-2.0 validation, a completed
  UMI-tools 1.1.2 log is a second explicit accepted shape. Its exact version,
  command/options, INFO record grammar, completion marker and count consistency
  are validated; only timestamps, host, pid/UUID and runtime provenance are
  ignored. The exact 1.x zero/zero empty-sample form (which has one additional
  blank line) is also accepted; other trailing prose or whitespace is not.
- **`CONFIG`** — `outdir` and optional `workDir` are validated non-empty
  run-specific paths before removal. The pinned pipeline's core parameter keys
  and their basic types must be present; that stable payload is the comparison
  contract.
- **`TRACE`** — the comparator requires the exact default trace columns emitted
  by pinned Nextflow 26.04.1. It ignores task IDs, hashes, native scheduler IDs,
  submission time, duration/realtime, CPU, memory and I/O measurements. It
  compares, with multiplicity and irrespective of row order, the derived
  process, full task name, tag, final status and exit code. `COMPLETED` and
  `CACHED` are both successful but remain distinct stable statuses. A failed or
  nonzero-exit task invalidates the trace even when both files contain it.
- **MultiQC data** — pinned MultiQC 1.14 writes a customer-published data
  directory. Its log includes wall-clock timestamps, Nextflow work paths,
  random temporary directories, and an external “new version available” check;
  those fields are validated then normalised, while module discovery, versions,
  levels and all other ordered messages remain stable. The JSON data object
  contains creation time and analysis/source paths; those provenance paths are
  removed or reduced to source basenames while general statistics, plot data,
  raw module data, title/version and other configuration remain exact.
  `multiqc_sources.txt` applies the same directory-only normalisation to its
  parsed source column. Citations and module data tables have no observed
  volatility and retain `BINARY_EXACT`.
- **`MTX`** — MatrixMarket entry ordering is not significant; the same matrix can
  be emitted in different row orders. Ordering differences are benign, so
  triplets are compared as a sorted multiset. **Value** differences at a given
  coordinate are real and fail outside envelope mode. Input is streamed into
  bounded-size external-sort runs with bounded merge fan-in; neither the
  decompressed archive nor all triplets/differences are retained in memory.
  Before comparison, each
  gzip member must decompress completely as UTF-8 and declare a MatrixMarket
  `matrix coordinate` object with `integer` or finite `real` values and
  `general` symmetry. Dimensions, declared nnz, entry arity, 1-based bounds,
  actual entry count and non-negative count values are all validated. Real
  values use exact decimal numeric canonicalisation (`1`, `1.0`, and `1e0`
  compare equal) without binary-float aliasing or underflow. Envelope sums use
  arbitrary-precision integers/exact decimals, so wrapping, rounding and
  floating overflow cannot turn a corrupt comparison into a pass. To keep that
  exact arithmetic resource-bounded on hostile input, integer tokens are
  limited to 1,000 digits, real tokens to 4,096 characters, and canonical real
  exponents to +/-10,000. Those limits are many orders of magnitude beyond a
  count-matrix producer value; exceeding one fails validation before summation.
- **`H5AD`** — byte identity is not meaningful for HDF5 containers. Both files
  must first open as AnnData, after which all customer-visible AnnData slots are
  compared. Count matrices are read in bounded row chunks and hashed in a
  storage-neutral canonical CSR form, so production-scale sparse matrices do
  not create a COO copy or one Python tuple per non-zero value. Empty,
  truncated, corrupt, non-AnnData, complex, negative, NaN, or infinite `X` data fail.
  Envelope row sums use Python integers or checked `math.fsum` to avoid silent
  numeric overflow. The pipeline deliberately emits zero-byte `*.empty.h5ad`
  sentinels for samples without a count matrix; only that suffix, with zero
  bytes on both sides, is accepted.
- **`BAM`** — compression and `@PG` invocation records are not alignment
  semantics. Each BAM must pass `samtools quickcheck` and be readable end to
  end. Within each `(reference, start)` tie, record order is immaterial and the
  comparator checks a multiset of the pipeline's stable molecule key:
  reference/start/strand, cell barcode, optional UMI (canonically empty for the
  current chemistry), and the complete `XT` gene assignment. XT absence is an
  explicit empty key for pre-annotation STAR and legitimately unassigned reads
  in the initial featureCounts BAM. It is required for every record in the
  high-confidence `*.mapped.sorted.filtered.annotated.bam` and downstream
  `*.dedup.bam`; symmetric tag loss in those outputs fails validation.
  QNAME source identity, SEQ, QUAL, MAPQ, CIGAR, mate fields, and
  incidental optional tags are excluded because parallel deduplication may
  select a different PCR duplicate to represent the same molecule. Comparison
  of coordinate-sorted files is bounded by the largest coordinate-tie group.
  The published STAR BAM is `SO:unsorted`, so it is first passed through
  samtools' 64-MiB external merge sort; this keeps Python memory bounded while
  using temporary disk proportional to that BAM.
  This exclusion is an intentional residual risk: a change confined to those
  representative-read fields will pass. Consumers who depend directly on
  QNAME/SEQ/QUAL/MAPQ/CIGAR/mate or incidental tags must add a separate
  comparison for their use case.
- **`HTML`** — Class D presentation artefacts contain volatile Plotly element
  IDs and run provenance, so their content is not an equivalence contract. They
  must nevertheless decode as UTF-8, parse as balanced HTML, and contain parsed
  elements and scripts expected for their published filename. Raw strings in
  HTML/JavaScript comments, JavaScript strings or text cannot satisfy a
  landmark, and error/failure headings are rejected. The consolidated report
  requires its title, main/overview/cross-sample landmarks and cross-sample
  table. MultiQC requires a report title, strictly parsed pinned-1.14 JSON in
  `script#mqc_config`, and
  `h1#page_title` elements, and at least one `mqc-module-section-*` section.
  Plotly outputs require a parsed graph div and a real `Plotly.newPlot` call
  targeting that div. Cell Caller deliberately emits empty plot sentinels for
  empty/degenerate samples; only `*_counts_pdf_with_threshold.html`,
  `*_barnyard_plot.html`, and `*_pdf_with_cutoff.html` (including the species
  variants) pass, and only when both sides are empty. Other empty HTML is a
  `DIFFER` and fails the run.

  The three `pipeline_info` HTML files have contracts taken from actual
  Nextflow 26.04.1 output. `execution_report.html` requires the report title,
  navigation, task table, and parsed `window.data.trace`/`summary`, rejects
  unsuccessful tasks, and compares the stable task and process-summary
  multisets.
  `execution_timeline.html` requires the timeline heading/elements and valid
  process timing records, ignores absolute timings, and compares process label
  plus cached state with multiplicity. `pipeline_dag.html` requires a parsed
  Mermaid flowchart and module initialiser and compares the flowchart topology
  and labels after excluding only the theme preamble.

### Published-tree audit

The configured `pipeline_info` outputs are all explicit: resolved configuration
is `CONFIG`, execution trace is `TRACE`, and report/timeline/DAG are `HTML` with
the Nextflow contracts above. Metrics and RSeQC text, dedup logs, raw/filtered
H5AD and tripartite matrices, STAR/featureCounts/dedup BAMs, consolidated/QC/
Cell Caller/MultiQC HTML, MultiQC log/JSON/sources, and the remaining stable
MultiQC data files likewise map to the table above. Any new HTML filename fails
closed until a contract is defined; any other new filename is included in the
file-set check and receives the fail-loud `BINARY_EXACT` fallback.
The five fixed `pipeline_info` files are required in both trees, so empty trees
or trees both missing a fixed manifest file cannot pass merely because their
file sets match. Dynamic sample outputs remain governed by symmetric file-set
equality; the comparator does not claim to infer their expected names without
the run's sample manifest.

For a 1.x-to-2.0 diagnostic, do not compare the two complete outdirs: 2.0
intentionally replaces report filenames and execution metadata. A curated
common-path inventory can be run with:

```bash
python3 compare_outputs.py <old_common> <new_common> \
  --allow-subset --envelope-max-flips N
```

Subset mode retains exact symmetric file-set checking and every per-class
contract; it disables only the five-file complete-run manifest requirement.
It is not expected to make the major-version comparison all green: documented
metric, RSeQC and dedup-count changes remain `DIFFER` and require review. The
strict default command is the future same-version reproducibility gate.

## Envelope mode (the count-matrix ambiguity)

The pipeline has one irreducible non-determinism. A tiny number of multimapped
reads are genuinely ambiguous between two genes, so a read can move between two
genes **for the same barcode** from one run to the next. This flips two
count-matrix entries but leaves the **per-barcode column total unchanged** — it
is *net-preserving*. This matrix invariant makes no promise about separately
derived metrics or RSeQC output: those classes remain exact, and their known
1.x-to-2.0 deltas correctly produce `DIFFER`.

`--envelope-max-flips N` turns on envelope mode for the `MTX` and `H5AD` classes
**only**. A count-matrix difference then passes with the distinct verdict
`ENVELOPE_OK` **iff all three** hold:

1. the MatrixMarket row/column shape and type metadata / H5AD `X` shape are
   identical (the validated MTX nnz may change when support changes);
2. **every per-column (per-barcode) sum is identical** between A and B; and
3. the number of differing `(row,col,value)` entries — the symmetric difference
   of the triplet multisets, counted as `max(only_in_a, only_in_b)` — is `<= N`.

Otherwise the verdict is `DIFFER` (fail). In particular, **a per-column-sum
change is always a `DIFFER`/FAIL regardless of `N`**: that is a real regression,
not the ambiguity envelope, because it changes a per-barcode total.

The **column-sum-preservation invariant** is the heart of the envelope: the
envelope only ever forgives differences that move counts *between genes within a
barcode*, never differences that change *how many counts a barcode has*.

Scope of relaxation: envelope mode relaxes **nothing else**. `TEXT_EXACT`
(metrics / multisample / RSeQC), `DEDUP_LOG`, `CONFIG`, `TRACE`, all MultiQC
classes, `GZ_TEXT_EXACT`
(barcodes / features), `BINARY_EXACT`, `BAM`, and `HTML` all retain their
documented strict contracts.

### Orientation: which axis is the barcode?

In `matrix.mtx.gz` the matrix is **genes-by-barcodes** (rows = genes, columns =
barcodes), so the per-column sum is the per-barcode total — summed over rows.

In the `.h5ad`, `adata.X` is the **transpose**: `obs`-by-`var` =
**barcodes-by-genes** (rows/`obs` = barcodes, columns/`var` = genes). To keep
the same "per-barcode total" semantics, the h5ad envelope sums `X` over `var`
(`axis=1`), i.e. per-`obs` sums. The two classes apply this invariant
independently. The comparator does not cross-link an H5AD to its companion
MatrixMarket/barcode/feature files, so it does not claim an H5AD-to-MTX
orientation or identity cross-check.

## Verdicts

`EQUAL`, `ENVELOPE_OK`, `DIFFER`, `PRESENT`, `MISSING` (the last only on
`FILE_SET_MISMATCH` rows). Only `DIFFER` in an exact class, or any
`FILE_SET_MISMATCH`, causes a nonzero exit. `EQUAL`, `ENVELOPE_OK` and `PRESENT`
pass. `ENVELOPE_OK` appears only when `--envelope-max-flips` is
set and is distinct from `EQUAL`: it records that a count-matrix difference was
*within* the characterised ambiguity envelope, not that the files were equal.

## Determinism probe

Comparing two runs of the same code on the same input should yield no failures
under these per-class contracts. That is semantic/output equivalence, not a
claim that the whole directory is byte-identical: HDF5/BAM containers,
presentation HTML and normalised Nextflow/MultiQC provenance deliberately have
other contracts. The raw and filtered tripartite count-matrix archives are the
stronger exception: from 2.0.0 their gzip headers and payload ordering are
deterministic, so identical inputs produce byte-identical `barcodes.tsv.gz`,
`features.tsv.gz`, and `matrix.mtx.gz` archives.

Note that the same upstream count differences surface in **both** the
`matrix.mtx.gz` and its `.h5ad` twin (the h5ad stores the same matrix,
barcodes-by-genes), so a genuine value difference is correctly reported in both.
