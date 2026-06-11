# Output-equivalence regression comparator

`compare_outputs.py` decides whether the published-output trees of two pipeline
runs are *output-equivalent*. It is foundational verification infrastructure for
a customer-facing pipeline: it fails loud and applies **no silent tolerance**.
A real divergence in any class that is meant to be exact drives a nonzero exit
code.

## Usage

```bash
python3 compare_outputs.py <dir_a> <dir_b> [--json OUT] [--max-diffs N]
```

- `<dir_a>`, `<dir_b>` — the two output directories to compare.
- `--json OUT` — also write the full per-file result set (with the run summary)
  as JSON to `OUT`.
- `--max-diffs N` — cap on the number of differing lines/triplets/keys reported
  per file (default 10). It only bounds the *report*, never the verdict.

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

## Verdicts

`EQUAL`, `DIFFER`, `SKIP`, `PRESENT`, `MISSING` (the last only on
`FILE_SET_MISMATCH` rows). Only `DIFFER` in an exact class, or any
`FILE_SET_MISMATCH`, causes a nonzero exit. `SKIP` and `PRESENT` never fail.

## Determinism probe

Comparing two runs of the same code on the same input should yield no failures.
Note that the same upstream count differences surface in **both** the
`matrix.mtx.gz` and its `.h5ad` twin (the h5ad stores the same matrix,
barcodes-by-genes), so a genuine value difference is correctly reported in both.
