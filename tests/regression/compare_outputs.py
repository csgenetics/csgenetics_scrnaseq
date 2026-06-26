#!/usr/bin/env python3
"""Output-equivalence regression comparator for the CS Genetics scRNA-seq pipeline.

Compares the published-output trees of two pipeline runs and decides, per file,
whether they are output-equivalent. Equivalence is class-specific: some files
must be byte-identical, some carry benign run-to-run volatility (umi_tools
timestamps, mtx row ordering, run-specific S3 paths) that must be normalised
away before comparison.

Design principles (see CLAUDE.md):
  - Fail loud. No silent tolerance: any real divergence in an exact class is a
    failure that drives a nonzero exit code.
  - The published-file SET is part of the contract; a missing/extra file fails.
  - Only normalise volatility that is provably run-specific and benign.

Usage:
    python3 compare_outputs.py <dir_a> <dir_b> [--json OUT] [--max-diffs N]
                               [--envelope-max-flips N]

Stdlib only, plus an OPTIONAL guarded `anndata` import for .h5ad files.

Envelope mode (--envelope-max-flips N): the pipeline has one irreducible
non-determinism. A tiny number of multimapped reads are genuinely ambiguous
between two genes, so a read can move between two genes for the SAME barcode
run-to-run. This flips two count-matrix entries but leaves the per-barcode
column total unchanged (net-preserving). Envelope mode PASSES a count-matrix
(MTX/H5AD) difference -- with the distinct verdict ENVELOPE_OK -- IFF the
per-column (per-barcode) sums are identical AND the number of differing entries
is within N. Every other class stays strict byte-exact; envelope mode never
relaxes them. A per-column-sum change is a real regression and always DIFFERs.
"""

import argparse
import gzip
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys

# ----------------------------------------------------------------------------
# Verdicts and class metadata
# ----------------------------------------------------------------------------

EQUAL = "EQUAL"
ENVELOPE_OK = "ENVELOPE_OK"  # passing: count-matrix diff inside the ambiguity envelope
DIFFER = "DIFFER"
SKIP = "SKIP"
PRESENT = "PRESENT"
MISSING = "MISSING"  # used only for FILE_SET_MISMATCH rows

# Verdicts that pass (never drive a nonzero exit).
PASSING_VERDICTS = {EQUAL, ENVELOPE_OK, SKIP, PRESENT}

# Classes whose DIFFER verdict is a real failure (must drive nonzero exit).
EXACT_CLASSES = {
    "TEXT_EXACT",
    "DEDUP_LOG",
    "CONFIG",
    "GZ_TEXT_EXACT",
    "MTX",
    "H5AD",
    "BINARY_EXACT",
    "BAM",
}


# ----------------------------------------------------------------------------
# Optional anndata
# ----------------------------------------------------------------------------

try:
    import anndata as _anndata  # noqa: F401
    _HAVE_ANNDATA = True
except Exception:
    _HAVE_ANNDATA = False

_HAVE_SAMTOOLS = shutil.which("samtools") is not None


# ----------------------------------------------------------------------------
# Classification
# ----------------------------------------------------------------------------

def classify(relpath):
    """Return the comparison class for a relative path by FIRST matching rule."""
    base = os.path.basename(relpath)
    # Normalise to forward slashes for the "/RSeQC/" containment test.
    norm = relpath.replace(os.sep, "/")

    if (
        relpath.endswith(".metrics.csv")
        or base == "multisample_out.csv"
        or ("/RSeQC/" in norm and relpath.endswith(".txt"))
    ):
        return "TEXT_EXACT"
    if relpath.endswith(".dedup.log"):
        return "DEDUP_LOG"
    if base == "resolved_configuration.txt":
        return "CONFIG"
    if relpath.endswith("barcodes.tsv.gz") or relpath.endswith("features.tsv.gz"):
        return "GZ_TEXT_EXACT"
    if relpath.endswith("matrix.mtx.gz"):
        return "MTX"
    if relpath.endswith(".h5ad"):
        return "H5AD"
    if relpath.endswith(".bam"):
        return "BAM"
    if relpath.endswith(".html"):
        return "HTML"
    return "BINARY_EXACT"


# ----------------------------------------------------------------------------
# Readers / helpers
# ----------------------------------------------------------------------------

def _read_bytes(path):
    with open(path, "rb") as fh:
        return fh.read()


def _read_text(path):
    with open(path, "r", encoding="utf-8", errors="surrogateescape") as fh:
        return fh.read()


def _gunzip_text(path):
    with gzip.open(path, "rt", encoding="utf-8", errors="surrogateescape") as fh:
        return fh.read()


def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _line_diffs(text_a, text_b, max_diffs):
    """Report up to max_diffs differing lines (1-based line numbers)."""
    a = text_a.splitlines()
    b = text_b.splitlines()
    diffs = []
    for i in range(max(len(a), len(b))):
        la = a[i] if i < len(a) else "<no line>"
        lb = b[i] if i < len(b) else "<no line>"
        if la != lb:
            diffs.append({"line": i + 1, "a": la, "b": lb})
            if len(diffs) >= max_diffs:
                break
    return diffs


# ----------------------------------------------------------------------------
# DEDUP_LOG
# ----------------------------------------------------------------------------

_DEDUP_PATTERNS = [
    re.compile(r"Input Reads:\s*(\d+)"),
    re.compile(r"Number of reads out:\s*(\d+)"),
    re.compile(r"Total number of positions deduplicated:\s*(\d+)"),
]


def _parse_dedup(path):
    """Extract the (input_reads, reads_out, positions_dedup) integer tuple.

    Handles two shapes:
      - a normal multi-line umi_tools log, and
      - the empty-sample branch, a single literal line containing escaped
        "\\n" sequences rather than real newlines.

    Each integer is required exactly once; anything else fails loud (returns an
    error marker) because a missing field means the log shape changed and the
    comparison can no longer be trusted.
    """
    raw = _read_text(path)
    # Normalise the literal-backslash-n empty-sample form into real newlines so
    # a single set of regexes covers both shapes.
    text = raw.replace("\\n", "\n")

    values = []
    for pat in _DEDUP_PATTERNS:
        matches = pat.findall(text)
        if len(matches) != 1:
            return None, (
                f"expected exactly one match for {pat.pattern!r}, "
                f"found {len(matches)}"
            )
        values.append(int(matches[0]))
    return tuple(values), None


def _compare_dedup(path_a, path_b, max_diffs):
    ta, ea = _parse_dedup(path_a)
    tb, eb = _parse_dedup(path_b)
    if ea or eb:
        detail = "; ".join(
            f"{side}: {err}"
            for side, err in (("a", ea), ("b", eb))
            if err
        )
        return DIFFER, {"parse_error": detail}
    if ta == tb:
        return EQUAL, {"tuple": list(ta)}
    return DIFFER, {"a": list(ta), "b": list(tb)}


# ----------------------------------------------------------------------------
# CONFIG (resolved_configuration.txt)
# ----------------------------------------------------------------------------

_CONFIG_DROP_KEYS = ("outdir", "workDir")


def _compare_config(path_a, path_b, max_diffs):
    try:
        da = json.loads(_read_text(path_a))
        db = json.loads(_read_text(path_b))
    except json.JSONDecodeError as exc:
        return DIFFER, {"parse_error": str(exc)}
    for key in _CONFIG_DROP_KEYS:
        da.pop(key, None)
        db.pop(key, None)
    if da == db:
        return EQUAL, {}
    diffs = []
    for key in sorted(set(da) | set(db)):
        va = da.get(key, "<absent>")
        vb = db.get(key, "<absent>")
        if va != vb:
            diffs.append({"key": key, "a": va, "b": vb})
            if len(diffs) >= max_diffs:
                break
    return DIFFER, {"diffs": diffs}


# ----------------------------------------------------------------------------
# Text / gz-text / binary
# ----------------------------------------------------------------------------

def _compare_text_exact(path_a, path_b, max_diffs):
    ba = _read_bytes(path_a)
    bb = _read_bytes(path_b)
    if ba == bb:
        return EQUAL, {}
    return DIFFER, {
        "diffs": _line_diffs(
            ba.decode("utf-8", "surrogateescape"),
            bb.decode("utf-8", "surrogateescape"),
            max_diffs,
        )
    }


def _compare_gz_text_exact(path_a, path_b, max_diffs):
    ta = _gunzip_text(path_a)
    tb = _gunzip_text(path_b)
    if ta == tb:
        return EQUAL, {}
    return DIFFER, {"diffs": _line_diffs(ta, tb, max_diffs)}


def _compare_binary_exact(path_a, path_b, max_diffs):
    if _sha256(path_a) == _sha256(path_b):
        return EQUAL, {}
    return DIFFER, {"sha256_a": _sha256(path_a), "sha256_b": _sha256(path_b)}


# ----------------------------------------------------------------------------
# MTX (Matrix Market)
# ----------------------------------------------------------------------------

def _parse_mtx(path):
    """Return (dims_header, sorted_triplet_lines).

    The MatrixMarket banner/comment lines (leading '%') are discarded; the first
    non-comment line is the dims header (compared exact); the remaining lines are
    the entry triplets whose ORDER is not significant.
    """
    text = _gunzip_text(path)
    dims = None
    triplets = []
    for line in text.splitlines():
        if line.startswith("%"):
            continue
        if dims is None:
            dims = line
            continue
        if line == "":
            continue
        triplets.append(line)
    triplets.sort()
    return dims, triplets


def _mtx_column_sums(triplets):
    """Map column index -> summed value over all rows, from MatrixMarket triplets.

    Each triplet is "row col value". In this pipeline's mtx the matrix is
    genes-by-barcodes (rows=genes, cols=barcodes), so the per-COLUMN sum is the
    per-barcode total -- the quantity the ambiguity envelope must preserve.
    """
    from collections import defaultdict

    sums = defaultdict(float)
    for line in triplets:
        row, col, val = line.split()
        sums[int(col)] += float(val)
    return dict(sums)


def _compare_mtx(path_a, path_b, max_diffs, envelope_max_flips=None):
    dims_a, trip_a = _parse_mtx(path_a)
    dims_b, trip_b = _parse_mtx(path_b)

    if dims_a != dims_b:
        return DIFFER, {"dims_a": dims_a, "dims_b": dims_b}

    if trip_a == trip_b:
        return EQUAL, {}

    # Sorted multiset symmetric difference. Walk both sorted lists in lockstep
    # so we can attribute each differing triplet to its side.
    from collections import Counter

    ca = Counter(trip_a)
    cb = Counter(trip_b)
    only_a = sorted((ca - cb).elements())
    only_b = sorted((cb - ca).elements())
    n_diff_entries = max(len(only_a), len(only_b))

    diffs = []
    for line in only_a:
        diffs.append({"side": "a", "triplet": line})
        if len(diffs) >= max_diffs:
            break
    if len(diffs) < max_diffs:
        for line in only_b:
            diffs.append({"side": "b", "triplet": line})
            if len(diffs) >= max_diffs:
                break

    detail = {
        "n_only_a": len(only_a),
        "n_only_b": len(only_b),
        "diffs": diffs,
    }

    if envelope_max_flips is None:
        return DIFFER, detail

    # Envelope mode: a difference passes IFF per-column (per-barcode) sums are
    # identical AND the differing-entry count is within the flip budget. A
    # column-sum change is a real regression -> always DIFFER regardless of N.
    sums_a = _mtx_column_sums(trip_a)
    sums_b = _mtx_column_sums(trip_b)
    diff_cols = sorted(
        c for c in (set(sums_a) | set(sums_b))
        if sums_a.get(c, 0.0) != sums_b.get(c, 0.0)
    )
    col_sums_preserved = not diff_cols
    detail.update({
        "n_diff_cols": len(diff_cols),
        "n_diff_entries": n_diff_entries,
        "col_sums_preserved": col_sums_preserved,
        "envelope_max_flips": envelope_max_flips,
    })

    if col_sums_preserved and n_diff_entries <= envelope_max_flips:
        return ENVELOPE_OK, detail
    return DIFFER, detail


# ----------------------------------------------------------------------------
# H5AD
# ----------------------------------------------------------------------------

def _h5ad_signature(path):
    """Build a comparable signature of an AnnData object's core content."""
    import numpy as np
    import scipy.sparse as sp

    adata = _anndata.read_h5ad(path)
    X = adata.X
    if sp.issparse(X):
        coo = X.tocoo()
        nnz = int(coo.nnz)
        entries = sorted(
            zip(coo.row.tolist(), coo.col.tolist(), coo.data.tolist())
        )
    else:
        arr = np.asarray(X)
        nz = np.nonzero(arr)
        nnz = int(len(nz[0]))
        entries = sorted(
            (int(i), int(j), arr[i, j])
            for i, j in zip(nz[0].tolist(), nz[1].tolist())
        )
    return {
        "shape": tuple(int(s) for s in X.shape),
        "nnz": nnz,
        "entries": entries,
        "var_names": list(adata.var_names),
        "obs_names": list(adata.obs_names),
    }


def _compare_h5ad(path_a, path_b, max_diffs, envelope_max_flips=None):
    if not _HAVE_ANNDATA:
        return SKIP, {"warn": "anndata not importable; .h5ad comparison skipped"}
    sig_a = _h5ad_signature(path_a)
    sig_b = _h5ad_signature(path_b)

    mismatches = {}
    if sig_a["shape"] != sig_b["shape"]:
        mismatches["shape"] = {"a": sig_a["shape"], "b": sig_b["shape"]}
    if sig_a["nnz"] != sig_b["nnz"]:
        mismatches["nnz"] = {"a": sig_a["nnz"], "b": sig_b["nnz"]}
    if sig_a["entries"] != sig_b["entries"]:
        mismatches["X_entries"] = "X (i,j,value) entries differ"
    # var/obs names: set AND order both matter.
    if sig_a["var_names"] != sig_b["var_names"]:
        mismatches["var_names"] = (
            "differ (set or order)"
            if set(sig_a["var_names"]) == set(sig_b["var_names"])
            else "differ (set)"
        )
    if sig_a["obs_names"] != sig_b["obs_names"]:
        mismatches["obs_names"] = (
            "differ (set or order)"
            if set(sig_a["obs_names"]) == set(sig_b["obs_names"])
            else "differ (set)"
        )

    if not mismatches:
        return EQUAL, {"shape": list(sig_a["shape"]), "nnz": sig_a["nnz"]}

    if envelope_max_flips is None:
        return DIFFER, {"mismatches": mismatches}

    # Envelope mode. adata.X is obs(barcodes)-by-var(genes) -- the TRANSPOSE of
    # the mtx (genes-by-barcodes). So the per-barcode total == the per-OBS sum
    # (sum over var, axis=1); that is the "column sum" / per-barcode invariant
    # the ambiguity envelope must preserve. The envelope only applies when the
    # difference is confined to X entries: any shape / name-axis mismatch is a
    # structural regression that always DIFFERs.
    structural = {k: v for k, v in mismatches.items() if k != "X_entries"}
    if structural:
        return DIFFER, {"mismatches": mismatches}

    ok, env_detail = _h5ad_envelope_check(
        path_a, path_b, sig_a, sig_b, envelope_max_flips
    )
    detail = {"mismatches": mismatches}
    detail.update(env_detail)
    return (ENVELOPE_OK if ok else DIFFER), detail


def _h5ad_envelope_check(path_a, path_b, sig_a, sig_b, envelope_max_flips):
    """Decide whether the X difference is inside the ambiguity envelope.

    Returns (passes, detail). Passes IFF per-barcode (per-obs) sums are identical
    AND the differing-entry count is within the flip budget.
    """
    import numpy as np

    # Differing X entries: symmetric difference of the (i,j,value) multisets,
    # counted as max(only_in_a, only_in_b). Entries are already sorted tuples.
    from collections import Counter

    ca = Counter(sig_a["entries"])
    cb = Counter(sig_b["entries"])
    n_only_a = sum((ca - cb).values())
    n_only_b = sum((cb - ca).values())
    n_diff_entries = max(n_only_a, n_only_b)

    # Per-obs (per-barcode) sums = sum over var (axis=1).
    adata_a = _anndata.read_h5ad(path_a)
    adata_b = _anndata.read_h5ad(path_b)
    sums_a = np.asarray(adata_a.X.sum(axis=1)).ravel()
    sums_b = np.asarray(adata_b.X.sum(axis=1)).ravel()
    col_sums_preserved = (
        sums_a.shape == sums_b.shape and np.array_equal(sums_a, sums_b)
    )
    n_diff_cols = (
        int(np.count_nonzero(sums_a != sums_b))
        if sums_a.shape == sums_b.shape
        else max(len(sums_a), len(sums_b))
    )

    detail = {
        "n_diff_cols": n_diff_cols,
        "n_diff_entries": n_diff_entries,
        "col_sums_preserved": col_sums_preserved,
        "envelope_max_flips": envelope_max_flips,
    }
    passes = col_sums_preserved and n_diff_entries <= envelope_max_flips
    return passes, detail


# ----------------------------------------------------------------------------
# BAM
# ----------------------------------------------------------------------------

def _samtools(args):
    return subprocess.run(
        ["samtools"] + args,
        check=True,
        capture_output=True,
        text=True,
    ).stdout


def _compare_bam(path_a, path_b, max_diffs):
    if not _HAVE_SAMTOOLS:
        return SKIP, {"warn": "samtools not on PATH; .bam comparison skipped"}
    flag_a = _samtools(["flagstat", path_a])
    flag_b = _samtools(["flagstat", path_b])
    count_a = _samtools(["view", "-c", path_a]).strip()
    count_b = _samtools(["view", "-c", path_b]).strip()

    if flag_a == flag_b and count_a == count_b:
        return EQUAL, {"count": count_a}
    detail = {}
    if flag_a != flag_b:
        detail["flagstat_diffs"] = _line_diffs(flag_a, flag_b, max_diffs)
    if count_a != count_b:
        detail["count"] = {"a": count_a, "b": count_b}
    return DIFFER, detail


# ----------------------------------------------------------------------------
# HTML
# ----------------------------------------------------------------------------

def _compare_html(path_a, path_b, max_diffs):
    size_a = os.path.getsize(path_a)
    size_b = os.path.getsize(path_b)
    # Class D presentation: both must exist and be non-empty; content not
    # compared. existence is guaranteed by pairing; verify size>0.
    if size_a > 0 and size_b > 0:
        return PRESENT, {"size_a": size_a, "size_b": size_b}
    # An empty published HTML is a real defect even though content is not
    # compared.
    return DIFFER, {"size_a": size_a, "size_b": size_b, "reason": "empty HTML"}


# ----------------------------------------------------------------------------
# Dispatch
# ----------------------------------------------------------------------------

_COMPARATORS = {
    "TEXT_EXACT": _compare_text_exact,
    "DEDUP_LOG": _compare_dedup,
    "CONFIG": _compare_config,
    "GZ_TEXT_EXACT": _compare_gz_text_exact,
    "MTX": _compare_mtx,
    "H5AD": _compare_h5ad,
    "BAM": _compare_bam,
    "HTML": _compare_html,
    "BINARY_EXACT": _compare_binary_exact,
}


# Classes whose comparator accepts the envelope flip budget. Every other class
# stays strict byte-exact; the envelope must NOT relax them.
_ENVELOPE_CLASSES = {"MTX", "H5AD"}


def compare_file(cls, path_a, path_b, max_diffs, envelope_max_flips=None):
    if cls in _ENVELOPE_CLASSES:
        return _COMPARATORS[cls](path_a, path_b, max_diffs, envelope_max_flips)
    return _COMPARATORS[cls](path_a, path_b, max_diffs)


# ----------------------------------------------------------------------------
# File listing
# ----------------------------------------------------------------------------

def list_relfiles(root):
    out = set()
    for dirpath, _dirnames, filenames in os.walk(root):
        for name in filenames:
            full = os.path.join(dirpath, name)
            out.add(os.path.relpath(full, root))
    return out


# ----------------------------------------------------------------------------
# Driver
# ----------------------------------------------------------------------------

def run(dir_a, dir_b, max_diffs, envelope_max_flips=None):
    files_a = list_relfiles(dir_a)
    files_b = list_relfiles(dir_b)

    only_a = sorted(files_a - files_b)
    only_b = sorted(files_b - files_a)
    paired = sorted(files_a & files_b)

    results = []
    file_set_mismatch = bool(only_a or only_b)

    for rel in only_a:
        results.append({
            "path": rel, "class": "FILE_SET_MISMATCH",
            "verdict": MISSING, "detail": {"present_in": "a"},
        })
    for rel in only_b:
        results.append({
            "path": rel, "class": "FILE_SET_MISMATCH",
            "verdict": MISSING, "detail": {"present_in": "b"},
        })

    for rel in paired:
        cls = classify(rel)
        verdict, detail = compare_file(
            cls, os.path.join(dir_a, rel), os.path.join(dir_b, rel),
            max_diffs, envelope_max_flips,
        )
        results.append({
            "path": rel, "class": cls, "verdict": verdict, "detail": detail,
        })

    results.sort(key=lambda r: r["path"])

    # Failure = any FILE_SET_MISMATCH or any DIFFER in an exact class. A passing
    # verdict (EQUAL/ENVELOPE_OK/SKIP/PRESENT) never fails; ENVELOPE_OK is the
    # envelope mode's passing verdict for an in-envelope count-matrix diff.
    failed = file_set_mismatch
    for r in results:
        if r["verdict"] == DIFFER and r["class"] in EXACT_CLASSES:
            failed = True

    summary = {
        "dir_a": os.path.abspath(dir_a),
        "dir_b": os.path.abspath(dir_b),
        "anndata_available": _HAVE_ANNDATA,
        "samtools_available": _HAVE_SAMTOOLS,
        "envelope_max_flips": envelope_max_flips,
        "file_set_mismatch": file_set_mismatch,
        "counts": _count_verdicts(results),
        "n_files_a": len(files_a),
        "n_files_b": len(files_b),
        "failed": failed,
    }
    return results, summary, failed


def _count_verdicts(results):
    counts = {}
    for r in results:
        counts[r["verdict"]] = counts.get(r["verdict"], 0) + 1
    return counts


# ----------------------------------------------------------------------------
# Reporting
# ----------------------------------------------------------------------------

def _format_envelope(detail):
    """One-line envelope summary, shared by MTX and H5AD."""
    return (
        f"      envelope: {detail['n_diff_entries']} differing entry(ies), "
        f"{detail['n_diff_cols']} differing column-sum(s), "
        f"column sums preserved={detail['col_sums_preserved']}, "
        f"max-flips={detail['envelope_max_flips']}"
    )


def _format_detail(r, max_diffs):
    cls = r["class"]
    detail = r["detail"]
    lines = []
    if r["verdict"] == ENVELOPE_OK and "col_sums_preserved" in detail:
        lines.append(_format_envelope(detail))
        return lines
    if cls in ("MTX", "H5AD") and r["verdict"] == DIFFER and "col_sums_preserved" in detail:
        lines.append(_format_envelope(detail))
    if cls == "MTX" and r["verdict"] == DIFFER:
        if "dims_a" in detail:
            lines.append(f"      dims differ: a={detail['dims_a']} b={detail['dims_b']}")
        if "n_only_a" in detail:
            lines.append(
                f"      {detail['n_only_a']} triplet(s) only in A, "
                f"{detail.get('n_only_b', '?')} only in B"
            )
        for d in detail.get("diffs", []):
            lines.append(f"      [{d['side']}] {d['triplet']}")
    elif cls == "DEDUP_LOG" and r["verdict"] == DIFFER:
        if "parse_error" in detail:
            lines.append(f"      parse error: {detail['parse_error']}")
        else:
            lines.append(f"      a={detail['a']} b={detail['b']}")
    elif r["verdict"] == DIFFER and "diffs" in detail:
        for d in detail["diffs"]:
            if "line" in d:
                lines.append(f"      line {d['line']}: a={d['a']!r} b={d['b']!r}")
            elif "key" in d:
                lines.append(f"      key {d['key']!r}: a={d['a']!r} b={d['b']!r}")
    elif cls in ("H5AD", "CONFIG", "BAM") and r["verdict"] == DIFFER:
        lines.append(f"      {json.dumps(detail)}")
    elif r["verdict"] in (SKIP, PRESENT) and detail.get("warn"):
        lines.append(f"      {detail['warn']}")
    return lines


def print_report(results, summary, max_diffs):
    print(f"Comparing:\n  A = {summary['dir_a']}\n  B = {summary['dir_b']}")
    emf = summary.get("envelope_max_flips")
    mode = "strict" if emf is None else f"envelope (max-flips={emf})"
    print(
        f"  anndata={'yes' if summary['anndata_available'] else 'no'}  "
        f"samtools={'yes' if summary['samtools_available'] else 'no'}  "
        f"mode={mode}"
    )
    print("-" * 78)
    pathw = max((len(r["path"]) for r in results), default=4)
    pathw = min(max(pathw, 4), 70)
    for r in results:
        print(f"  {r['verdict']:<8} {r['class']:<17} {r['path']}")
        for line in _format_detail(r, max_diffs):
            print(line)
    print("-" * 78)
    counts = summary["counts"]
    counts_str = "  ".join(f"{k}={v}" for k, v in sorted(counts.items()))
    print(f"SUMMARY: {counts_str}")
    if summary["file_set_mismatch"]:
        print("FILE_SET_MISMATCH: file sets differ between A and B")
    print(f"RESULT: {'FAIL' if summary['failed'] else 'PASS'}")


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Output-equivalence regression comparator for the "
                    "CS Genetics scRNA-seq pipeline."
    )
    parser.add_argument("dir_a")
    parser.add_argument("dir_b")
    parser.add_argument("--json", dest="json_out", default=None,
                        help="write the full result set as JSON to this path")
    parser.add_argument("--max-diffs", type=int, default=10,
                        help="max differing lines/triplets to report per file")
    parser.add_argument("--envelope-max-flips", type=int, default=None,
                        help="enable count-matrix envelope mode: PASS (ENVELOPE_OK) "
                             "an MTX/H5AD diff iff per-barcode column sums are "
                             "identical and the differing-entry count is <= this "
                             "budget. Default (unset) keeps strict byte-exact mode. "
                             "All other classes stay strict regardless.")
    args = parser.parse_args(argv)

    for d in (args.dir_a, args.dir_b):
        if not os.path.isdir(d):
            parser.error(f"not a directory: {d}")

    if args.envelope_max_flips is not None and args.envelope_max_flips < 0:
        parser.error("--envelope-max-flips must be >= 0")

    results, summary, failed = run(
        args.dir_a, args.dir_b, args.max_diffs, args.envelope_max_flips
    )
    print_report(results, summary, args.max_diffs)

    if args.json_out:
        with open(args.json_out, "w", encoding="utf-8") as fh:
            json.dump({"summary": summary, "results": results}, fh, indent=2)
        print(f"JSON written to {args.json_out}")

    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
