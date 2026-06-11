# Output-equivalence regression comparator

`compare_outputs.py` decides whether the published-output trees of two pipeline
runs are *output-equivalent*. It is foundational verification infrastructure for
a customer-facing pipeline: it fails loud and applies **no silent tolerance**.
A real divergence in any class that is meant to be exact drives a nonzero exit
code.

## Usage

```bash
python3 compare_outputs.py <dir_a> <dir_b> [--json OUT] [--max-diffs N] \
                           [--envelope-max-flips N]
```

- `<dir_a>`, `<dir_b>` — the two output directories to compare.
- `--json OUT` — also write the full per-file result set (with the run summary)
  as JSON to `OUT`.
- `--max-diffs N` — cap on the number of differing lines/triplets/keys reported
  per file (default 10). It only bounds the *report*, never the verdict.
- `--envelope-max-flips N` — enable **envelope mode** for the count matrices
  (`MTX`, `H5AD`). Unset (default) = strict byte-exact mode. See
  [Envelope mode](#envelope-mode-the-count-matrix-ambiguity).

Exit code is `0` when equivalent, `1` when any failure is detected.

Dependencies: Python standard library only. `anndata` is an **optional** import;
if it is not importable, `.h5ad` files are reported `SKIP` (with a warning)
rather than failing. `samtools` on `PATH` is likewise optional for `.bam`.

## Contract: the published-file set

The *set* of relative paths is part of the contract. If any relative path exists
under one directory but not the other, that is a `FILE_SET_MISMATCH` and the run
fails. Extra or missing published files are never tolerated.

## Per-class rules

Each paired file is classified by the **first** matching rule below, then
compared by that class's logic. Only provably run-specific, benign volatility is
normalised away; everything else is compared exactly.

| Class | Matches | Comparison | DIFFER = failure? |
| --- | --- | --- | --- |
| `TEXT_EXACT` | `*.metrics.csv`, `multisample_out.csv`, `*.txt` under a `/RSeQC/` path | exact byte compare; reports up to `--max-diffs` differing lines | yes |
| `DEDUP_LOG` | `*.dedup.log` | extract the three integers (`Input Reads`, `Number of reads out`, `Total number of positions deduplicated`) and compare only those tuples | yes |
| `CONFIG` | `resolved_configuration.txt` | parse pretty-printed JSON, drop `outdir` and `workDir`, compare the remaining dict | yes |
| `GZ_TEXT_EXACT` | `*barcodes.tsv.gz`, `*features.tsv.gz` | gunzip, exact compare | yes |
| `MTX` | `*matrix.mtx.gz` | gunzip; dims header compared exact; entry triplets compared as a **sorted multiset** | yes |
| `H5AD` | `*.h5ad` | (anndata) compare `X` (shape, nnz, sorted `(i,j,value)` entries), `var_names` and `obs_names` (set **and** order) | yes (SKIP if anndata absent) |
| `BAM` | `*.bam` | (samtools) compare `flagstat` text and `view -c` count; no byte compare | yes (SKIP if samtools absent) |
| `HTML` | `*.html` | presentation only: assert both exist and are non-empty; content not compared | no (`PRESENT`) |
| `BINARY_EXACT` | anything else | sha256 compare | yes |

### Why these normalisations

- **`DEDUP_LOG`** — umi_tools logs carry volatile timestamps, pid, an instance
  UUID, and CPU-time lines that change every run. Only the read-count integers
  are deterministic output. The empty-sample branch writes a single literal line
  with escaped `\n` rather than real newlines; both shapes parse to the same
  three integers.
- **`CONFIG`** — `outdir` and `workDir` are run-specific S3 paths; the rest of
  the resolved params is contract.
- **`MTX`** — MatrixMarket entry ordering is not significant; the same matrix can
  be emitted in different row orders. Ordering differences are benign, so
  triplets are compared as a sorted multiset. **Value** differences at a given
  coordinate are real and fail. This is the key rule.
- **`HTML`** — Class D presentation artefacts; their content is not part of the
  numeric contract, so only presence and non-emptiness are checked.

## Envelope mode (the count-matrix ambiguity)

The pipeline has one irreducible non-determinism. A tiny number of multimapped
reads are genuinely ambiguous between two genes, so a read can move between two
genes **for the same barcode** from one run to the next. This flips two
count-matrix entries but leaves the **per-barcode column total unchanged** — it
is *net-preserving*. The customer metrics CSVs (`*.metrics.csv`,
`multisample_out.csv`, RSeQC) stay byte-stable; only the raw count matrix
(`matrix.mtx.gz` and its `.h5ad` twin) differs.

`--envelope-max-flips N` turns on envelope mode for the `MTX` and `H5AD` classes
**only**. A count-matrix difference then passes with the distinct verdict
`ENVELOPE_OK` **iff all three** hold:

1. the dims header / `X` shape is identical;
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
(metrics / multisample / RSeQC), `DEDUP_LOG`, `CONFIG`, `GZ_TEXT_EXACT`
(barcodes / features), `BINARY_EXACT` and `BAM` all stay strict byte-exact.

### Orientation: which axis is the barcode?

In `matrix.mtx.gz` the matrix is **genes-by-barcodes** (rows = genes, columns =
barcodes), so the per-column sum is the per-barcode total — summed over rows.

In the `.h5ad`, `adata.X` is the **transpose**: `obs`-by-`var` =
**barcodes-by-genes** (rows/`obs` = barcodes, columns/`var` = genes). To keep
the same "per-barcode total" semantics, the h5ad envelope sums `X` over `var`
(`axis=1`), i.e. per-`obs` sums. This is verified against the mtx dims so the
two classes test the identical invariant.

## Verdicts

`EQUAL`, `ENVELOPE_OK`, `DIFFER`, `SKIP`, `PRESENT`, `MISSING` (the last only on
`FILE_SET_MISMATCH` rows). Only `DIFFER` in an exact class, or any
`FILE_SET_MISMATCH`, causes a nonzero exit. `EQUAL`, `ENVELOPE_OK`, `SKIP` and
`PRESENT` all pass. `ENVELOPE_OK` appears only when `--envelope-max-flips` is
set and is distinct from `EQUAL`: it records that a count-matrix difference was
*within* the characterised ambiguity envelope, not that the files were equal.

## Determinism probe

Comparing two runs of the same code on the same input should yield no failures.
Note that the same upstream count differences surface in **both** the
`matrix.mtx.gz` and its `.h5ad` twin (the h5ad stores the same matrix,
barcodes-by-genes), so a genuine value difference is correctly reported in both.
