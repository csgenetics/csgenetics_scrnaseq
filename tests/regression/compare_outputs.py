#!/usr/bin/env python3
"""Output-equivalence regression comparator for the CS Genetics scRNA-seq pipeline.

Compares the published-output trees of two pipeline runs and decides, per file,
whether they are output-equivalent. Equivalence is class-specific: some files
must be byte-identical, some carry benign run-to-run volatility (umi_tools
timestamps, MatrixMarket entry ordering, run-specific S3/work paths, Nextflow
runtime measurements and MultiQC provenance) that must be normalised away
before comparison.

Design principles (see CLAUDE.md):
  - Fail loud. No silent tolerance: any real divergence in an exact class is a
    failure that drives a nonzero exit code.
  - The published-file SET is part of the contract; a missing/extra file fails.
  - Only normalise volatility that is provably run-specific and benign.

Usage:
    python3 compare_outputs.py <dir_a> <dir_b> [--json OUT] [--max-diffs N]
                               [--envelope-max-flips N] [--allow-subset]

Stdlib only unless the compared tree contains H5AD or BAM files. Those classes
require `anndata` or `samtools`, respectively, and fail closed if the required
reader is unavailable.

Envelope mode (--envelope-max-flips N): the pipeline has one irreducible
non-determinism. A tiny number of multimapped reads are genuinely ambiguous
between two genes, so a read can move between two genes for the SAME barcode
run-to-run. This flips two count-matrix entries but leaves the per-barcode
column total unchanged (net-preserving). Envelope mode PASSES a count-matrix
(MTX/H5AD) difference -- with the distinct verdict ENVELOPE_OK -- IFF the
per-column (per-barcode) sums are identical AND the number of differing entries
is within N. Every other class keeps its documented strict contract; envelope
mode never relaxes them. A per-column-sum change is a real regression and always
DIFFERs.
"""

import argparse
import csv
import decimal
import gzip
import hashlib
import heapq
import html.parser
import itertools
import json
import math
import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile

# ----------------------------------------------------------------------------
# Verdicts and class metadata
# ----------------------------------------------------------------------------

EQUAL = "EQUAL"
ENVELOPE_OK = "ENVELOPE_OK"  # passing: count-matrix diff inside the ambiguity envelope
DIFFER = "DIFFER"
PRESENT = "PRESENT"
MISSING = "MISSING"  # used only for FILE_SET_MISMATCH rows

# Classes whose DIFFER verdict is a real failure (must drive nonzero exit).
EXACT_CLASSES = {
    "TEXT_EXACT",
    "DEDUP_LOG",
    "CONFIG",
    "TRACE",
    "MULTIQC_LOG",
    "MULTIQC_JSON",
    "MULTIQC_SOURCES",
    "GZ_TEXT_EXACT",
    "MTX",
    "H5AD",
    "BINARY_EXACT",
    "BAM",
    "HTML",
}

# These files are emitted for every completed pipeline run, independently of
# sample count. Pairwise set equality alone would otherwise let two empty trees
# or two trees missing the same fixed run-manifest artifact pass.
_REQUIRED_PUBLISHED_FILES = frozenset({
    "pipeline_info/resolved_configuration.txt",
    "pipeline_info/execution_trace.txt",
    "pipeline_info/execution_report.html",
    "pipeline_info/execution_timeline.html",
    "pipeline_info/pipeline_dag.html",
})


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
        or ((norm.startswith("RSeQC/") or "/RSeQC/" in norm)
            and relpath.endswith(".txt"))
    ):
        return "TEXT_EXACT"
    if relpath.endswith(".dedup.log"):
        return "DEDUP_LOG"
    if base == "resolved_configuration.txt":
        return "CONFIG"
    if base == "execution_trace.txt":
        return "TRACE"
    if base == "multiqc.log":
        return "MULTIQC_LOG"
    if base == "multiqc_data.json" or base.endswith(".multiqc.data.json"):
        return "MULTIQC_JSON"
    if base == "multiqc_sources.txt":
        return "MULTIQC_SOURCES"
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
    with open(path, "r", encoding="utf-8", errors="strict") as fh:
        return fh.read()


def _gunzip_text(path):
    with gzip.open(path, "rt", encoding="utf-8", errors="strict") as fh:
        return fh.read()


def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _line_diffs(text_a, text_b, max_diffs):
    """Report up to max_diffs differing lines (1-based line numbers)."""
    if max_diffs == 0:
        return []
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


def _reject_json_constant(value):
    raise ValueError(f"non-finite JSON constant {value}")


def _json_object_without_duplicates(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate JSON object key {key!r}")
        result[key] = value
    return result


def _validate_finite_json(value, location="$"):
    if isinstance(value, float) and not math.isfinite(value):
        raise ValueError(f"non-finite JSON number at {location}")
    if isinstance(value, dict):
        for key, item in value.items():
            _validate_finite_json(item, f"{location}.{key}")
    elif isinstance(value, list):
        for index, item in enumerate(value):
            _validate_finite_json(item, f"{location}[{index}]")


def _load_json_strict(text):
    value = json.loads(
        text,
        parse_constant=_reject_json_constant,
        object_pairs_hook=_json_object_without_duplicates,
    )
    _validate_finite_json(value)
    return value


# ----------------------------------------------------------------------------
# DEDUP_LOG
# ----------------------------------------------------------------------------

_DEDUP_LOG_RE = re.compile(
    r"\AINFO Reads: Input Reads: (\d+)\n"
    r"INFO Number of reads out: (\d+)\n\Z"
)
_DEDUP_LEGACY_EMPTY_LOG = (
    "INFO Reads: Input Reads: 0\n"
    "INFO Number of reads out: 0\n\n"
)
_DEDUP_LEGACY_INFO_RE = re.compile(
    r"\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2},\d{3} INFO (?P<message>.+)"
)
_DEDUP_LEGACY_STABLE_OPTIONS = {
    "assigned_tag": "None",
    "cell_tag": "None",
    "cell_tag_delim": "None",
    "cell_tag_split": "-",
    "chimeric_pairs": "use",
    "chrom": "None",
    "compresslevel": "6",
    "detection_method": "None",
    "filter_umi": "None",
    "gene_tag": "None",
    "gene_transcript_map": "None",
    "get_umi_method": "read_id",
    "ignore_tlen": "False",
    "ignore_umi": "False",
    "in_sam": "True",
    "log2stderr": "False",
    "loglevel": "1",
    "mapping_quality": "0",
    "method": "directional",
    "no_sort_output": "False",
    "out_sam": "False",
    "output_unmapped": "False",
    "paired": "False",
    "per_cell": "True",
    "per_contig": "False",
    "per_gene": "False",
    "random_seed": "None",
    "read_length": "False",
    "short_help": "None",
    "skip_regex": "^(__|Unassigned)",
    "soft_clip_threshold": "4",
    "spliced": "False",
    "stats": "False",
    "stderr": "<_io.TextIOWrapper name='<stderr>' mode='w' encoding='utf-8'>",
    "stdout": "<_io.TextIOWrapper name='<stdout>' mode='w' encoding='utf-8'>",
    "subset": "None",
    "threshold": "1",
    "timeit_file": "None",
    "timeit_header": "None",
    "timeit_name": "all",
    "tmpdir": "None",
    "umi_sep": "_",
    "umi_tag": "RX",
    "umi_tag_delim": "None",
    "umi_tag_split": "None",
    "umi_whitelist": "None",
    "umi_whitelist_paired": "None",
    "unmapped_reads": "discard",
    "unpaired_reads": "use",
    "whole_contig": "False",
}
_DEDUP_LEGACY_DYNAMIC_OPTIONS = {"stdin", "stdlog"}


def _single_legacy_dedup_metric(messages, pattern, label):
    matches = [re.fullmatch(pattern, message) for message in messages]
    values = [match.group(1) for match in matches if match]
    if len(values) != 1:
        raise ValueError(
            f"legacy umi_tools log must contain exactly one {label}, found {len(values)}"
        )
    return int(values[0])


def _parse_legacy_dedup(raw):
    """Validate the completed UMI-tools 1.1.2 log published by pipeline 1.x."""
    lines = raw.splitlines()
    if len(lines) < 12 or lines[0] != "# UMI-tools version: 1.1.2":
        raise ValueError("not a UMI-tools 1.1.2 legacy log")
    generated = re.fullmatch(
        r"# output generated by dedup --per-cell --in-sam "
        r"-I (?P<input>\S+) --log=(?P<log>\S+)",
        lines[1],
    )
    if not generated:
        raise ValueError("legacy umi_tools log has an unexpected generated-command header")
    input_name = generated.group("input")
    log_name = generated.group("log")
    if not re.fullmatch(r"# job started at .+ on \S+ -- \S+", lines[2]):
        raise ValueError("legacy umi_tools log lacks its job-start record")
    if not re.fullmatch(r"# pid: \d+, system: .+", lines[3]):
        raise ValueError("legacy umi_tools log lacks its pid/system record")
    if not re.fullmatch(r"# job finished in \d+ seconds at .+ -- .+ -- \S+", lines[-1]):
        raise ValueError("legacy umi_tools log lacks its completion record")

    options = {}
    messages = []
    for line_number, line in enumerate(lines[4:-1], start=5):
        option = re.fullmatch(r"# ([a-z][a-z0-9_]*)\s+: (.*)", line)
        if option:
            key, value = option.groups()
            if key in options:
                raise ValueError(
                    f"legacy umi_tools log line {line_number}: duplicate option {key}"
                )
            options[key] = value
            continue
        info = _DEDUP_LEGACY_INFO_RE.fullmatch(line)
        if info:
            messages.append(info.group("message"))
            continue
        raise ValueError(
            f"legacy umi_tools log line {line_number}: unexpected record shape"
        )
    expected_keys = set(_DEDUP_LEGACY_STABLE_OPTIONS) | _DEDUP_LEGACY_DYNAMIC_OPTIONS
    if set(options) != expected_keys:
        raise ValueError(
            "legacy umi_tools option set differs: "
            f"missing={sorted(expected_keys - set(options))}, "
            f"unexpected={sorted(set(options) - expected_keys)}"
        )
    for key, expected in _DEDUP_LEGACY_STABLE_OPTIONS.items():
        if options[key] != expected:
            raise ValueError(
                f"legacy umi_tools option {key!r} must be {expected!r}, "
                f"found {options[key]!r}"
            )
    expected_stdin = (
        f"<_io.TextIOWrapper name='{input_name}' mode='r' encoding='UTF-8'>"
    )
    expected_stdlog = (
        f"<_io.TextIOWrapper name='{log_name}' mode='a' encoding='UTF-8'>"
    )
    if options["stdin"] != expected_stdin or options["stdlog"] != expected_stdlog:
        raise ValueError("legacy umi_tools stream options do not match command paths")
    commands = [
        message for message in messages
        if message.startswith("command: ")
    ]
    expected_command = (
        f"command: dedup --per-cell --in-sam -I {input_name} --log={log_name}"
    )
    if commands != [expected_command]:
        raise ValueError("legacy umi_tools log has an invalid command record")
    allowed_messages = (
        r"command: dedup --per-cell --in-sam -I \S+ --log=\S+",
        r"Written out \d+ reads",
        r"Parsed \d+ input reads",
        r"Reads: Input Reads: \d+",
        r"Number of reads out: \d+",
        r"Total number of positions deduplicated: \d+",
        r"Mean number of unique UMIs per position: \d+(?:\.\d+)?",
        r"Max\. number of unique UMIs per position: \d+",
    )
    for message in messages:
        if not any(re.fullmatch(pattern, message) for pattern in allowed_messages):
            raise ValueError(f"legacy umi_tools log has unexpected INFO record {message!r}")

    input_reads = _single_legacy_dedup_metric(
        messages, r"Reads: Input Reads: (\d+)", "Input Reads field"
    )
    reads_out = _single_legacy_dedup_metric(
        messages, r"Number of reads out: (\d+)", "reads-out field"
    )
    positions = _single_legacy_dedup_metric(
        messages,
        r"Total number of positions deduplicated: (\d+)",
        "positions-deduplicated field",
    )
    if positions != reads_out:
        raise ValueError(
            "legacy umi_tools positions-deduplicated count must equal reads out "
            f"({positions} != {reads_out})"
        )
    parsed_progress = [
        int(match.group(1))
        for message in messages
        if (match := re.fullmatch(r"Parsed (\d+) input reads", message))
    ]
    written_progress = [
        int(match.group(1))
        for message in messages
        if (match := re.fullmatch(r"Written out (\d+) reads", message))
    ]
    expected_parsed = list(range(1_000_000, input_reads + 1, 1_000_000))
    expected_written = list(range(100_000, reads_out + 1, 100_000))
    if parsed_progress != expected_parsed:
        raise ValueError(
            "legacy umi_tools parsed-read progress is inconsistent: "
            f"expected {expected_parsed}, found {parsed_progress}"
        )
    if written_progress != expected_written:
        raise ValueError(
            "legacy umi_tools written-read progress is inconsistent: "
            f"expected {expected_written}, found {written_progress}"
        )
    means = [
        match.group(1)
        for message in messages
        if (match := re.fullmatch(
            r"Mean number of unique UMIs per position: (\d+(?:\.\d+)?)",
            message,
        ))
    ]
    maxima = [
        match.group(1)
        for message in messages
        if (match := re.fullmatch(
            r"Max\. number of unique UMIs per position: (\d+)", message
        ))
    ]
    if means != ["1.00"] or maxima != ["1"]:
        raise ValueError(
            "legacy umi_tools UMI summary must be exactly mean=1.00, max=1"
        )
    return input_reads, reads_out


def _parse_dedup(path):
    """Extract the current published ``(input_reads, reads_out)`` tuple.

    Pipeline 2.0 producer branches write the same exact two-line summary after
    volatile per-contig logs have been consolidated. For documented 1.x->2.0
    comparisons, the completed UMI-tools 1.1.2 log is a second explicit shape;
    only its known timestamp/host/pid/runtime provenance is ignored.
    """
    raw = _read_text(path)
    match = _DEDUP_LOG_RE.fullmatch(raw)
    if match:
        values = [int(value) for value in match.groups()]
    elif raw == _DEDUP_LEGACY_EMPTY_LOG:
        values = [0, 0]
    else:
        try:
            values = list(_parse_legacy_dedup(raw))
        except ValueError as exc:
            return None, (
                "expected the exact 2.0 two-line summary or a validated completed "
                f"UMI-tools 1.1.2 log: {exc}"
            )
    input_reads, reads_out = values
    if reads_out > input_reads:
        return None, (
            "Number of reads out cannot exceed Input Reads "
            f"({reads_out} > {input_reads})"
        )
    return tuple(values), None


def _compare_dedup(path_a, path_b, max_diffs):
    parsed = {}
    read_errors = {}
    for side, path in (("a", path_a), ("b", path_b)):
        try:
            parsed[side] = _parse_dedup(path)
        except (OSError, UnicodeError) as exc:
            read_errors[side] = f"{type(exc).__name__}: {exc}"
    if read_errors:
        return DIFFER, {"validation_errors": read_errors}
    ta, ea = parsed["a"]
    tb, eb = parsed["b"]
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
_CONFIG_REQUIRED_STABLE_KEYS = frozenset({
    "input_csv",
    "barcode_list_path",
    "barcode_correction_list_path",
    "barcode_pattern",
    "minimum_count_threshold",
    "sss_nmer",
    "awsregion",
    "awsqueue",
    "max_memory",
    "max_cpus",
    "max_time",
})


def _normalise_config(data):
    """Validate the published Nextflow parameter map before dropping paths."""
    if not isinstance(data, dict) or not data:
        raise ValueError(
            "resolved configuration must be a non-empty JSON object"
        )
    if "outdir" not in data:
        raise ValueError("resolved configuration lacks required outdir")
    for key in _CONFIG_DROP_KEYS:
        if key not in data:
            continue
        if not isinstance(data[key], str) or not data[key].strip():
            raise ValueError(
                f"resolved configuration {key} must be a non-empty string"
            )
    stable = {key: value for key, value in data.items()
              if key not in _CONFIG_DROP_KEYS}
    missing = sorted(_CONFIG_REQUIRED_STABLE_KEYS - stable.keys())
    if missing:
        raise ValueError(
            f"resolved configuration lacks stable pipeline parameters: {missing}"
        )
    for key in ("barcode_list_path", "barcode_correction_list_path"):
        if not isinstance(stable[key], str) or not stable[key].strip():
            raise ValueError(f"resolved configuration {key} must be non-empty")
    if not isinstance(stable["barcode_pattern"], str) \
            or not re.fullmatch(r"C+", stable["barcode_pattern"]):
        raise ValueError(
            "resolved configuration barcode_pattern must contain only C bases"
        )
    threshold = stable["minimum_count_threshold"]
    if isinstance(threshold, bool) or not isinstance(threshold, (int, float)) \
            or (isinstance(threshold, float) and not math.isfinite(threshold)) \
            or threshold < 0:
        raise ValueError(
            "resolved configuration minimum_count_threshold must be "
            "a finite non-negative number"
        )
    for key in ("sss_nmer", "max_cpus"):
        if isinstance(stable[key], bool) or not isinstance(stable[key], int) \
                or stable[key] <= 0:
            raise ValueError(
                f"resolved configuration {key} must be a positive integer"
            )
    for key in ("max_memory", "max_time"):
        if not isinstance(stable[key], dict) or not stable[key]:
            raise ValueError(
                f"resolved configuration {key} must be a non-empty object"
            )
    for key in ("input_csv", "awsregion", "awsqueue"):
        if stable[key] is not None and (
            not isinstance(stable[key], str) or not stable[key].strip()
        ):
            raise ValueError(
                f"resolved configuration {key} must be null or a non-empty string"
            )
    return stable


def _compare_config(path_a, path_b, max_diffs):
    try:
        da = _load_json_strict(_read_text(path_a))
        db = _load_json_strict(_read_text(path_b))
    except (OSError, UnicodeError, ValueError) as exc:
        return DIFFER, {"parse_error": str(exc)}
    try:
        da = _normalise_config(da)
        db = _normalise_config(db)
    except ValueError as exc:
        return DIFFER, {"parse_error": str(exc)}
    if da == db:
        return EQUAL, {}
    diffs = []
    if max_diffs == 0:
        return DIFFER, {"diffs": diffs}
    for key in sorted(set(da) | set(db)):
        va = da.get(key, "<absent>")
        vb = db.get(key, "<absent>")
        if va != vb:
            diffs.append({"key": key, "a": va, "b": vb})
            if len(diffs) >= max_diffs:
                break
    return DIFFER, {"diffs": diffs}


# ----------------------------------------------------------------------------
# TRACE (Nextflow 26.04.1 execution_trace.txt)
# ----------------------------------------------------------------------------

_TRACE_COLUMNS = (
    "task_id",
    "hash",
    "native_id",
    "name",
    "status",
    "exit",
    "submit",
    "duration",
    "realtime",
    "%cpu",
    "peak_rss",
    "peak_vmem",
    "rchar",
    "wchar",
)

_TRACE_VOLATILE_COLUMNS = {
    "task_id",
    "hash",
    "native_id",
    "submit",
    "duration",
    "realtime",
    "%cpu",
    "peak_rss",
    "peak_vmem",
    "rchar",
    "wchar",
}

_TRACE_FINAL_STATUSES = {"COMPLETED", "CACHED", "FAILED", "ABORTED"}
_TRACE_SUCCESS_STATUSES = {"COMPLETED", "CACHED"}
_TRACE_TIMESTAMP_RE = re.compile(r"\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}(?:\.\d+)?")
_TRACE_DURATION_RE = re.compile(
    r"(?:-|(?:\d+(?:\.\d+)?(?:ms|s|m|h|d))(?: \d+(?:\.\d+)?(?:ms|s|m|h|d))*)"
)
_TRACE_CPU_RE = re.compile(r"(?:-|\d+(?:\.\d+)?%)")
_TRACE_SIZE_RE = re.compile(r"(?:-|0|\d+(?:\.\d+)? (?:B|KB|MB|GB|TB|PB))")
_TRACE_VOLATILE_PATTERNS = {
    "submit": _TRACE_TIMESTAMP_RE,
    "duration": _TRACE_DURATION_RE,
    "realtime": _TRACE_DURATION_RE,
    "%cpu": _TRACE_CPU_RE,
    "peak_rss": _TRACE_SIZE_RE,
    "peak_vmem": _TRACE_SIZE_RE,
    "rchar": _TRACE_SIZE_RE,
    "wchar": _TRACE_SIZE_RE,
}


def _split_trace_task_name(name):
    """Return the process and optional tag encoded in Nextflow's task name."""
    if name.endswith(")") and " (" in name:
        process, _separator, suffix = name.partition(" (")
        if process:
            return process, suffix[:-1]
    return name, ""


def _parse_trace(path):
    """Parse the exact default trace schema emitted by pinned Nextflow 26.04.1."""
    try:
        with open(path, "r", encoding="utf-8", errors="strict", newline="") as fh:
            rows = list(csv.reader(fh, delimiter="\t", strict=True))
    except (OSError, UnicodeError, csv.Error) as exc:
        raise ValueError(f"cannot read trace TSV: {exc}") from exc

    if not rows:
        raise ValueError("trace is empty")
    header = tuple(rows[0])
    if header != _TRACE_COLUMNS:
        raise ValueError(
            "unexpected trace header; expected "
            f"{list(_TRACE_COLUMNS)!r}, found {list(header)!r}"
        )
    if len(rows) == 1:
        raise ValueError("trace contains no task rows")

    records = []
    failures = []
    task_ids = set()
    for line_number, values in enumerate(rows[1:], start=2):
        if len(values) != len(_TRACE_COLUMNS):
            raise ValueError(
                f"line {line_number}: expected {len(_TRACE_COLUMNS)} columns, "
                f"found {len(values)}"
            )
        row = dict(zip(_TRACE_COLUMNS, values))
        if any(row[column] == "" for column in _TRACE_COLUMNS):
            raise ValueError(f"line {line_number}: empty trace field")
        try:
            task_id = int(row["task_id"])
            exit_code = int(row["exit"])
        except ValueError as exc:
            raise ValueError(
                f"line {line_number}: task_id and exit must be integers"
            ) from exc
        if task_id < 1 or exit_code < 0:
            raise ValueError(
                f"line {line_number}: task_id must be positive and exit non-negative"
            )
        if task_id in task_ids:
            raise ValueError(f"line {line_number}: duplicate task_id {task_id}")
        task_ids.add(task_id)
        if not re.fullmatch(r"[0-9a-f]{2}/[0-9a-f]{6}", row["hash"]):
            raise ValueError(f"line {line_number}: invalid Nextflow task hash")
        if row["status"] not in _TRACE_FINAL_STATUSES:
            raise ValueError(
                f"line {line_number}: unsupported final status {row['status']!r}"
            )
        if not row["name"].strip():
            raise ValueError(f"line {line_number}: empty task name")
        for column, pattern in _TRACE_VOLATILE_PATTERNS.items():
            if not pattern.fullmatch(row[column]):
                raise ValueError(
                    f"line {line_number}: invalid {column} value {row[column]!r}"
                )

        process, tag = _split_trace_task_name(row["name"])
        record = (process, row["name"], tag, row["status"], exit_code)
        records.append(record)
        if row["status"] not in _TRACE_SUCCESS_STATUSES or exit_code != 0:
            failures.append(record)
    return records, failures


def _trace_record_json(record, multiplicity=1):
    process, task, tag, status, exit_code = record
    return {
        "process": process,
        "task": task,
        "tag": tag,
        "status": status,
        "exit": exit_code,
        "multiplicity": multiplicity,
    }


def _compare_trace(path_a, path_b, max_diffs):
    from collections import Counter

    parsed = {}
    validation_errors = {}
    for side, path in (("a", path_a), ("b", path_b)):
        try:
            parsed[side] = _parse_trace(path)
        except ValueError as exc:
            validation_errors[side] = str(exc)
    if validation_errors:
        return DIFFER, {"validation_errors": validation_errors}

    records_a, failures_a = parsed["a"]
    records_b, failures_b = parsed["b"]
    counter_a = Counter(records_a)
    counter_b = Counter(records_b)
    detail = {
        "task_count_a": len(records_a),
        "task_count_b": len(records_b),
        "normalised_fields": sorted(_TRACE_VOLATILE_COLUMNS),
    }
    if failures_a or failures_b:
        detail["execution_failures"] = {
            "a": [_trace_record_json(record) for record in failures_a[:max_diffs]],
            "b": [_trace_record_json(record) for record in failures_b[:max_diffs]],
        }
    if counter_a != counter_b:
        detail["only_a"] = [
            _trace_record_json(record, count)
            for record, count in sorted((counter_a - counter_b).items())[:max_diffs]
        ]
        detail["only_b"] = [
            _trace_record_json(record, count)
            for record, count in sorted((counter_b - counter_a).items())[:max_diffs]
        ]
    if failures_a or failures_b or counter_a != counter_b:
        return DIFFER, detail
    return EQUAL, {"task_count": len(records_a)}


# ----------------------------------------------------------------------------
# MultiQC 1.14 published data/provenance
# ----------------------------------------------------------------------------

_MULTIQC_LOG_RE = re.compile(
    r"^\[(?P<timestamp>\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2},\d{3})\] "
    r"(?P<logger>\S(?:.*?\S)?)\s+\[(?P<level>[A-Z]+)\s*\]\s+"
    r"(?P<message>.*)$"
)
_MULTIQC_ERROR_LEVELS = {"ERROR", "CRITICAL"}


def _normalise_multiqc_command(command):
    try:
        tokens = shlex.split(command)
    except ValueError as exc:
        raise ValueError(f"malformed MultiQC command: {exc}") from exc
    if not tokens:
        raise ValueError("empty MultiQC command")
    tokens[0] = os.path.basename(tokens[0])
    return shlex.join(tokens)


def _normalise_multiqc_log_message(message):
    if re.fullmatch(r"MultiQC Version v\S+ now available!", message):
        # External update-check result, unrelated to the pinned pipeline run.
        return None
    if message.startswith("Command used: "):
        command = _normalise_multiqc_command(message.removeprefix("Command used: "))
        return f"Command used: {command}"
    for prefix, replacement in (
        ("Working dir : ", "Working dir : <RUN_DIR>"),
        ("Search path : ", "Search path : <RUN_DIR>"),
        (
            "Using temporary directory for creating report: ",
            "Using temporary directory for creating report: <TEMP_DIR>",
        ),
    ):
        if message.startswith(prefix):
            value = message.removeprefix(prefix)
            if not value.strip() or "\x00" in value:
                raise ValueError(
                    f"MultiQC {prefix.strip()} path must be non-empty"
                )
            return replacement
    for label in ("Report", "Data"):
        match = re.fullmatch(rf"{label}\s+:\s+(.+)", message)
        if match:
            value = match.group(1)
            if not value.strip() or "\x00" in value:
                raise ValueError(f"MultiQC {label} path must be non-empty")
            filename = os.path.basename(value.replace("\\", "/"))
            if not filename or filename in {".", ".."}:
                raise ValueError(f"MultiQC {label} path lacks a filename")
            return f"{label}: {filename}"
    move = re.fullmatch(r"Moving data file from '([^']+)' to '([^']+)'", message)
    if move:
        source, target = move.groups()
        if any(not value.strip() or "\x00" in value for value in (source, target)):
            raise ValueError("MultiQC moved-data paths must be non-empty")
        destination = os.path.basename(target.replace("\\", "/"))
        if not destination or destination in {".", ".."}:
            raise ValueError("MultiQC moved-data destination lacks a filename")
        return f"Moving data file from <TEMP_DIR>/multiqc_data to <RUN_DIR>/{destination}"
    return message


def _parse_multiqc_log(path):
    try:
        lines = _read_text(path).splitlines()
    except (OSError, UnicodeError) as exc:
        raise ValueError(f"cannot read MultiQC log: {exc}") from exc
    if not lines:
        raise ValueError("MultiQC log is empty")

    records = []
    has_version = False
    has_completion = False
    for line_number, line in enumerate(lines, start=1):
        match = _MULTIQC_LOG_RE.fullmatch(line)
        if not match:
            raise ValueError(f"line {line_number}: malformed MultiQC log record")
        logger = match.group("logger").strip()
        level = match.group("level")
        message = match.group("message")
        if level in _MULTIQC_ERROR_LEVELS:
            raise ValueError(
                f"line {line_number}: MultiQC logged {level}: {message}"
            )
        if message == "This is MultiQC v1.14":
            has_version = True
        if message == "MultiQC complete":
            has_completion = True
        message = _normalise_multiqc_log_message(message)
        if message is not None:
            records.append((logger, level, message))
    if not has_version:
        raise ValueError("MultiQC log lacks the pinned v1.14 version record")
    if not has_completion:
        raise ValueError("MultiQC log lacks its completion record")
    return records


def _compare_multiqc_log(path_a, path_b, max_diffs):
    parsed = {}
    validation_errors = {}
    for side, path in (("a", path_a), ("b", path_b)):
        try:
            parsed[side] = _parse_multiqc_log(path)
        except ValueError as exc:
            validation_errors[side] = str(exc)
    if validation_errors:
        return DIFFER, {"validation_errors": validation_errors}
    if parsed["a"] == parsed["b"]:
        return EQUAL, {"record_count": len(parsed["a"])}
    return DIFFER, {
        "diffs": _line_diffs(
            "\n".join(map(repr, parsed["a"])),
            "\n".join(map(repr, parsed["b"])),
            max_diffs,
        )
    }


_MULTIQC_JSON_REQUIRED = {
    "report_data_sources",
    "report_general_stats_data",
    "report_general_stats_headers",
    "report_multiqc_command",
    "report_plot_data",
    "report_saved_raw_data",
    "config_analysis_dir_abs",
    "config_analysis_dir",
    "config_creation_date",
    "config_title",
    "config_version",
}
_MULTIQC_JSON_VOLATILE = {
    "config_analysis_dir_abs",
    "config_analysis_dir",
    "config_creation_date",
}


def _normalise_multiqc_source_paths(value):
    if isinstance(value, dict):
        if not value:
            raise ValueError(
                "MultiQC report_data_sources containers must not be empty"
            )
        if any(not key.strip() for key in value):
            raise ValueError(
                "MultiQC report_data_sources keys must not be blank"
            )
        return {
            key: _normalise_multiqc_source_paths(item)
            for key, item in value.items()
        }
    if isinstance(value, list):
        if not value:
            raise ValueError(
                "MultiQC report_data_sources containers must not be empty"
            )
        return [_normalise_multiqc_source_paths(item) for item in value]
    if not isinstance(value, str):
        raise ValueError("MultiQC report_data_sources leaves must be paths")
    if not value.strip() or "\x00" in value:
        raise ValueError(
            "MultiQC report_data_sources paths must be non-empty strings"
        )
    basename = os.path.basename(value.replace("\\", "/"))
    if not basename or basename in {".", ".."}:
        raise ValueError(
            "MultiQC report_data_sources paths must have a source basename"
        )
    return basename


def _parse_multiqc_json(path):
    try:
        data = _load_json_strict(_read_text(path))
    except (OSError, UnicodeError, ValueError) as exc:
        raise ValueError(f"cannot parse MultiQC JSON: {exc}") from exc
    if not isinstance(data, dict):
        raise ValueError("MultiQC JSON must be an object")
    missing = sorted(_MULTIQC_JSON_REQUIRED - data.keys())
    if missing:
        raise ValueError(f"MultiQC JSON lacks required fields: {missing}")
    if not isinstance(data["report_data_sources"], dict):
        raise ValueError("MultiQC report_data_sources must be an object")
    if not data["report_data_sources"]:
        raise ValueError("MultiQC report_data_sources must not be empty")
    if not isinstance(data["report_general_stats_data"], list):
        raise ValueError("MultiQC report_general_stats_data must be a list")
    if not isinstance(data["report_general_stats_headers"], (list, dict)):
        raise ValueError("MultiQC report_general_stats_headers has invalid type")
    for key in ("report_plot_data", "report_saved_raw_data"):
        if not isinstance(data[key], dict):
            raise ValueError(f"MultiQC {key} must be an object")
    for key in ("report_multiqc_command", "config_creation_date", "config_title", "config_version"):
        if not isinstance(data[key], str) or not data[key].strip():
            raise ValueError(f"MultiQC {key} must be a non-empty string")
    if data["config_version"] != "1.14":
        raise ValueError(
            f"unexpected MultiQC version {data['config_version']!r}; expected 1.14"
        )
    for key in ("config_analysis_dir_abs", "config_analysis_dir"):
        paths = data[key]
        if not isinstance(paths, list) or not paths:
            raise ValueError(f"MultiQC {key} must be a non-empty path list")
        for path in paths:
            if not isinstance(path, str) or not path.strip() or "\x00" in path:
                raise ValueError(
                    f"MultiQC {key} entries must be non-empty strings"
                )

    data["report_data_sources"] = _normalise_multiqc_source_paths(
        data["report_data_sources"]
    )
    data["report_multiqc_command"] = _normalise_multiqc_command(
        data["report_multiqc_command"]
    )
    for key in _MULTIQC_JSON_VOLATILE:
        data.pop(key)
    return data


def _compare_multiqc_json(path_a, path_b, max_diffs):
    parsed = {}
    validation_errors = {}
    for side, path in (("a", path_a), ("b", path_b)):
        try:
            parsed[side] = _parse_multiqc_json(path)
        except ValueError as exc:
            validation_errors[side] = str(exc)
    if validation_errors:
        return DIFFER, {"validation_errors": validation_errors}
    if parsed["a"] == parsed["b"]:
        return EQUAL, {}
    differing_keys = [
        key for key in sorted(set(parsed["a"]) | set(parsed["b"]))
        if parsed["a"].get(key) != parsed["b"].get(key)
    ]
    return DIFFER, {"differing_keys": differing_keys[:max_diffs]}


_MULTIQC_SOURCES_HEADER = ("Module", "Section", "Sample Name", "Source")


def _parse_multiqc_sources(path):
    try:
        with open(path, "r", encoding="utf-8", errors="strict", newline="") as fh:
            rows = list(csv.reader(fh, delimiter="\t", strict=True))
    except (OSError, UnicodeError, csv.Error) as exc:
        raise ValueError(f"cannot parse MultiQC sources TSV: {exc}") from exc
    if not rows or tuple(rows[0]) != _MULTIQC_SOURCES_HEADER:
        raise ValueError("unexpected MultiQC sources TSV header")
    if len(rows) == 1:
        raise ValueError("MultiQC sources TSV contains no source rows")
    records = []
    for line_number, row in enumerate(rows[1:], start=2):
        if len(row) != len(_MULTIQC_SOURCES_HEADER) or any(
            not value.strip() or "\x00" in value for value in row
        ):
            raise ValueError(f"line {line_number}: invalid MultiQC source record")
        module, section, sample, source = row
        records.append(
            (module, section, sample, _normalise_multiqc_source_paths(source))
        )
    return records


def _compare_multiqc_sources(path_a, path_b, max_diffs):
    from collections import Counter

    parsed = {}
    validation_errors = {}
    for side, path in (("a", path_a), ("b", path_b)):
        try:
            parsed[side] = Counter(_parse_multiqc_sources(path))
        except ValueError as exc:
            validation_errors[side] = str(exc)
    if validation_errors:
        return DIFFER, {"validation_errors": validation_errors}
    if parsed["a"] == parsed["b"]:
        return EQUAL, {"source_count": sum(parsed["a"].values())}
    return DIFFER, {
        "only_a": [list(record) + [count]
                   for record, count in (parsed["a"] - parsed["b"]).items()][:max_diffs],
        "only_b": [list(record) + [count]
                   for record, count in (parsed["b"] - parsed["a"]).items()][:max_diffs],
    }


# ----------------------------------------------------------------------------
# Text / gz-text / binary
# ----------------------------------------------------------------------------

def _compare_text_exact(path_a, path_b, max_diffs):
    ba = _read_bytes(path_a)
    bb = _read_bytes(path_b)
    texts = {}
    validation_errors = {}
    for side, content in (("a", ba), ("b", bb)):
        try:
            texts[side] = content.decode("utf-8", "strict")
        except UnicodeDecodeError as exc:
            validation_errors[side] = f"invalid UTF-8: {exc}"
    if validation_errors:
        return DIFFER, {"validation_errors": validation_errors}
    if ba == bb:
        return EQUAL, {}
    return DIFFER, {
        "diffs": _line_diffs(texts["a"], texts["b"], max_diffs)
    }


def _compare_gz_text_exact(path_a, path_b, max_diffs):
    texts = {}
    validation_errors = {}
    for side, path in (("a", path_a), ("b", path_b)):
        try:
            texts[side] = _gunzip_text(path)
        except (OSError, EOFError, UnicodeError) as exc:
            validation_errors[side] = f"{type(exc).__name__}: {exc}"
    if validation_errors:
        return DIFFER, {"validation_errors": validation_errors}
    ta = texts["a"]
    tb = texts["b"]
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

_MTX_SORT_CHUNK_BYTES = 8 << 20
_MTX_MERGE_FAN_IN = 32
_MTX_MAX_INTEGER_DIGITS = 1_000
_MTX_MAX_DECIMAL_TOKEN_CHARS = 4_096
_MTX_MAX_DECIMAL_EXPONENT = 10_000
_MTX_REAL_RE = re.compile(
    r"[+-]?(?:(?:\d+(?:\.\d*)?)|(?:\.\d+))(?:[eE][+-]?\d+)?"
)


def _canonical_decimal_token(token, line_number):
    """Return an exact, context-free coefficient/exponent representation."""
    if len(token) > _MTX_MAX_DECIMAL_TOKEN_CHARS:
        raise ValueError(
            f"line {line_number}: real value exceeds the supported "
            f"{_MTX_MAX_DECIMAL_TOKEN_CHARS}-character safety limit"
        )
    if token.lower().lstrip("+-") in {"nan", "snan", "inf", "infinity"}:
        raise ValueError(f"line {line_number}: non-finite matrix value")
    if not _MTX_REAL_RE.fullmatch(token):
        raise ValueError(f"line {line_number}: invalid real value")
    try:
        value = decimal.Decimal(token)
    except decimal.InvalidOperation as exc:
        raise ValueError(f"line {line_number}: invalid real value") from exc
    if not value.is_finite():
        raise ValueError(f"line {line_number}: non-finite matrix value")
    if value < 0:
        raise ValueError(
            f"line {line_number}: count-matrix values must be non-negative"
        )
    if value.is_zero():
        return "0"
    sign, digits, exponent = value.as_tuple()
    if sign:
        raise ValueError(
            f"line {line_number}: count-matrix values must be non-negative"
        )
    digits = list(digits)
    while digits and digits[-1] == 0:
        digits.pop()
        exponent += 1
    if abs(exponent) > _MTX_MAX_DECIMAL_EXPONENT:
        raise ValueError(
            f"line {line_number}: real exponent exceeds the supported "
            f"+/-{_MTX_MAX_DECIMAL_EXPONENT} exact-sum safety limit"
        )
    coefficient = "".join(str(digit) for digit in digits)
    return f"{coefficient}e{exponent}"


def _parse_mtx_integer_token(token, line_number, label):
    """Parse one strict integer without relying on interpreter digit limits."""
    if not re.fullmatch(r"[+-]?\d+", token):
        raise ValueError(f"line {line_number}: {label} must be an integer")
    digit_count = len(token.lstrip("+-"))
    if digit_count > _MTX_MAX_INTEGER_DIGITS:
        raise ValueError(
            f"line {line_number}: {label} exceeds the supported "
            f"{_MTX_MAX_INTEGER_DIGITS}-digit safety limit"
        )
    try:
        return int(token)
    except ValueError as exc:
        raise ValueError(
            f"line {line_number}: invalid integer {label}"
        ) from exc


def _spill_mtx_run(records, directory, sequence):
    records.sort()
    path = os.path.join(directory, f"run-{sequence:08d}.txt")
    with open(path, "w", encoding="ascii", newline="") as fh:
        fh.writelines(records)
    return path


def _merge_mtx_run_group(paths, output):
    handles = [open(path, "r", encoding="ascii", newline="") for path in paths]
    try:
        with open(output, "w", encoding="ascii", newline="") as destination:
            destination.writelines(heapq.merge(*handles))
    finally:
        for handle in handles:
            handle.close()


def _collapse_mtx_runs(paths, directory, sequence):
    """Bound merge fan-in so very large archives cannot exhaust file handles."""
    paths = list(paths)
    generation = 0
    while len(paths) > _MTX_MERGE_FAN_IN:
        collapsed = []
        for offset in range(0, len(paths), _MTX_MERGE_FAN_IN):
            group = paths[offset:offset + _MTX_MERGE_FAN_IN]
            if len(group) == 1:
                collapsed.extend(group)
                continue
            output = os.path.join(
                directory, f"merge-{generation:04d}-{sequence:08d}.txt"
            )
            sequence += 1
            _merge_mtx_run_group(group, output)
            for path in group:
                os.unlink(path)
            collapsed.append(output)
        paths = collapsed
        generation += 1
    return paths


def _prepare_mtx(path, directory):
    """Validate an MTX archive and external-sort canonical records by column.

    Memory is bounded by one input line plus ``_MTX_SORT_CHUNK_BYTES`` of
    canonical records (and their fixed per-record Python overhead). Sorted runs
    spill to the comparator-owned temporary directory.
    """
    os.makedirs(directory, exist_ok=True)
    dimensions = None
    field = None
    object_type = storage = symmetry = None
    entry_count = 0
    records = []
    record_bytes = 0
    runs = []
    sequence = 0

    def spill():
        nonlocal records, record_bytes, sequence
        if records:
            runs.append(_spill_mtx_run(records, directory, sequence))
            sequence += 1
            records = []
            record_bytes = 0

    try:
        with gzip.open(
            path, "rt", encoding="utf-8", errors="strict", newline=""
        ) as fh:
            first = fh.readline()
            if first == "":
                raise ValueError("empty MatrixMarket file")
            banner = first.strip().split()
            if len(banner) != 5 or banner[0] != "%%MatrixMarket":
                raise ValueError("invalid MatrixMarket banner")
            object_type, storage, field, symmetry = (
                token.lower() for token in banner[1:]
            )
            if object_type != "matrix":
                raise ValueError(f"unsupported MatrixMarket object {object_type!r}")
            if storage != "coordinate":
                raise ValueError(f"unsupported MatrixMarket format {storage!r}")
            if field not in {"integer", "real"}:
                raise ValueError(f"unsupported MatrixMarket field {field!r}")
            if symmetry != "general":
                raise ValueError(f"unsupported MatrixMarket symmetry {symmetry!r}")

            for line_number, line in enumerate(fh, start=2):
                stripped = line.strip()
                if not stripped or stripped.startswith("%"):
                    continue
                tokens = stripped.split()
                if dimensions is None:
                    if len(tokens) != 3:
                        raise ValueError(
                            f"line {line_number}: dimensions must contain three integers"
                        )
                    n_rows, n_cols, declared_nnz = (
                        _parse_mtx_integer_token(
                            token, line_number, "MatrixMarket dimension"
                        )
                        for token in tokens
                    )
                    if n_rows < 0 or n_cols < 0 or declared_nnz < 0:
                        raise ValueError(
                            f"line {line_number}: dimensions and declared nnz "
                            "must be non-negative"
                        )
                    dimensions = (n_rows, n_cols, declared_nnz)
                    row_width = max(1, len(str(n_rows)))
                    column_width = max(1, len(str(n_cols)))
                    continue

                if len(tokens) != 3:
                    raise ValueError(
                        f"line {line_number}: coordinate entry must contain "
                        "row, column, value"
                    )
                row = _parse_mtx_integer_token(tokens[0], line_number, "row")
                column = _parse_mtx_integer_token(
                    tokens[1], line_number, "column"
                )
                n_rows, n_cols, _declared_nnz = dimensions
                if not (1 <= row <= n_rows) or not (1 <= column <= n_cols):
                    raise ValueError(
                        f"line {line_number}: coordinate ({row}, {column}) outside "
                        f"1-based bounds ({n_rows}, {n_cols})"
                    )

                if field == "integer":
                    numeric = _parse_mtx_integer_token(
                        tokens[2], line_number, "integer value"
                    )
                    if numeric < 0:
                        raise ValueError(
                            f"line {line_number}: count-matrix values must be "
                            "non-negative"
                        )
                    value_token = str(numeric)
                else:
                    value_token = _canonical_decimal_token(tokens[2], line_number)

                record = (
                    f"{column:0{column_width}d}\t{row:0{row_width}d}\t"
                    f"{value_token}\n"
                )
                records.append(record)
                record_bytes += len(record)
                entry_count += 1
                if record_bytes >= _MTX_SORT_CHUNK_BYTES:
                    spill()
    except (OSError, EOFError, UnicodeError) as exc:
        raise ValueError(
            f"cannot read gzip/UTF-8 MatrixMarket data: {exc}"
        ) from exc

    if dimensions is None:
        raise ValueError("missing MatrixMarket dimensions")
    if entry_count != dimensions[2]:
        raise ValueError(
            f"declared nnz is {dimensions[2]}, but parsed {entry_count} entries"
        )
    spill()
    runs = _collapse_mtx_runs(runs, directory, sequence)
    return {
        "object": object_type,
        "format": storage,
        "field": field,
        "symmetry": symmetry,
        "shape": dimensions[:2],
        "declared_nnz": dimensions[2],
        "runs": tuple(runs),
    }


def _iter_mtx_records(matrix):
    handles = [
        open(path, "r", encoding="ascii", newline="")
        for path in matrix["runs"]
    ]
    try:
        yield from heapq.merge(*handles)
    finally:
        for handle in handles:
            handle.close()


def _mtx_record_parts(record):
    column, row, value = record.rstrip("\n").split("\t")
    return int(row), int(column), value


def _format_mtx_record(record):
    row, column, value = _mtx_record_parts(record)
    return f"{row} {column} {value}"


def _compare_mtx_record_streams(matrix_a, matrix_b, max_diffs):
    iterator_a = iter(_iter_mtx_records(matrix_a))
    iterator_b = iter(_iter_mtx_records(matrix_b))
    record_a = next(iterator_a, None)
    record_b = next(iterator_b, None)
    n_only_a = n_only_b = 0
    diffs = []
    while record_a is not None or record_b is not None:
        if record_a is not None and record_a == record_b:
            record_a = next(iterator_a, None)
            record_b = next(iterator_b, None)
        elif record_b is None or (
            record_a is not None and record_a < record_b
        ):
            n_only_a += 1
            if len(diffs) < max_diffs:
                diffs.append({
                    "side": "a", "triplet": _format_mtx_record(record_a)
                })
            record_a = next(iterator_a, None)
        else:
            n_only_b += 1
            if len(diffs) < max_diffs:
                diffs.append({
                    "side": "b", "triplet": _format_mtx_record(record_b)
                })
            record_b = next(iterator_b, None)
    return n_only_a, n_only_b, diffs


def _add_exact_decimal(total, value_token):
    """Add a canonical non-negative decimal without context rounding."""
    if value_token == "0":
        return total
    coefficient_text, _separator, exponent_text = value_token.partition("e")
    coefficient = int(coefficient_text)
    exponent = int(exponent_text)
    total_coefficient, total_exponent = total
    if total_coefficient == 0:
        return coefficient, exponent
    common_exponent = min(total_exponent, exponent)
    if max(total_exponent, exponent) - common_exponent \
            > 2 * _MTX_MAX_DECIMAL_EXPONENT:
        raise ValueError("real-value exponent span exceeds exact-sum safety limit")
    total_coefficient = (
        total_coefficient * (10 ** (total_exponent - common_exponent))
        + coefficient * (10 ** (exponent - common_exponent))
    )
    total_exponent = common_exponent
    while total_coefficient and total_coefficient % 10 == 0:
        total_coefficient //= 10
        total_exponent += 1
    return total_coefficient, total_exponent


def _iter_mtx_column_sums(matrix):
    current_column = None
    total = 0 if matrix["field"] == "integer" else (0, 0)
    for record in _iter_mtx_records(matrix):
        _row, column, value_token = _mtx_record_parts(record)
        if current_column is not None and column != current_column:
            yield current_column, total
            total = 0 if matrix["field"] == "integer" else (0, 0)
        current_column = column
        if matrix["field"] == "integer":
            total += int(value_token)
        else:
            total = _add_exact_decimal(total, value_token)
    if current_column is not None:
        yield current_column, total


def _count_mtx_column_sum_differences(matrix_a, matrix_b):
    iterator_a = iter(_iter_mtx_column_sums(matrix_a))
    iterator_b = iter(_iter_mtx_column_sums(matrix_b))
    item_a = next(iterator_a, None)
    item_b = next(iterator_b, None)
    zero = 0 if matrix_a["field"] == "integer" else (0, 0)
    differences = 0
    while item_a is not None or item_b is not None:
        if item_a is not None and item_b is not None and item_a[0] == item_b[0]:
            differences += item_a[1] != item_b[1]
            item_a = next(iterator_a, None)
            item_b = next(iterator_b, None)
        elif item_b is None or (item_a is not None and item_a[0] < item_b[0]):
            differences += item_a[1] != zero
            item_a = next(iterator_a, None)
        else:
            differences += item_b[1] != zero
            item_b = next(iterator_b, None)
    return differences


def _compare_mtx(path_a, path_b, max_diffs, envelope_max_flips=None):
    with tempfile.TemporaryDirectory(prefix="compare-mtx-") as scratch:
        matrices = {}
        validation_errors = {}
        for side, path in (("a", path_a), ("b", path_b)):
            try:
                matrices[side] = _prepare_mtx(
                    path, os.path.join(scratch, side)
                )
            except ValueError as exc:
                validation_errors[side] = str(exc)
        if validation_errors:
            return DIFFER, {"validation_errors": validation_errors}

        matrix_a = matrices["a"]
        matrix_b = matrices["b"]
        # Declared nnz is validated against each actual stream, but envelope
        # mode may legitimately accept a support change with a different nnz.
        metadata_keys = ("object", "format", "field", "symmetry", "shape")
        metadata_a = {key: matrix_a[key] for key in metadata_keys}
        metadata_b = {key: matrix_b[key] for key in metadata_keys}
        if metadata_a != metadata_b:
            return DIFFER, {"metadata_a": metadata_a, "metadata_b": metadata_b}

        n_only_a, n_only_b, diffs = _compare_mtx_record_streams(
            matrix_a, matrix_b, max_diffs
        )
        if n_only_a == 0 and n_only_b == 0:
            return EQUAL, {}

        n_diff_entries = max(n_only_a, n_only_b)
        detail = {
            "n_only_a": n_only_a,
            "n_only_b": n_only_b,
            "diffs": diffs,
        }
        if envelope_max_flips is None:
            return DIFFER, detail

        try:
            n_diff_cols = _count_mtx_column_sum_differences(matrix_a, matrix_b)
        except ValueError as exc:
            detail["validation_errors"] = {"envelope": str(exc)}
            return DIFFER, detail
        col_sums_preserved = n_diff_cols == 0
        detail.update({
            "n_diff_cols": n_diff_cols,
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

def _normalise_scalar(value):
    """Return a type-aware, hashable representation of a scalar value."""
    import numpy as np
    import pandas as pd

    if isinstance(value, np.generic):
        value = value.item()
    if value is None or value is pd.NA:
        return ("missing",)
    if isinstance(value, float):
        if math.isnan(value):
            return ("float", "nan")
        if math.isinf(value):
            return ("float", "inf" if value > 0 else "-inf")
        return ("float", value)
    if isinstance(value, bool):
        return ("bool", value)
    if isinstance(value, int):
        return ("int", value)
    if isinstance(value, complex):
        return ("complex", value.real, value.imag)
    if isinstance(value, str):
        return ("str", value)
    if isinstance(value, bytes):
        return ("bytes", value.hex())

    # Timestamp, timedelta and other extension scalars have stable string
    # representations, but keep their type so unlike values cannot collide.
    return (type(value).__name__, str(value))


_MATRIX_CHUNK_TARGET_BYTES = 4 << 20


def _digest_scalar_sequence(values):
    """Hash a scalar sequence without retaining Python objects for every value."""
    digest = hashlib.sha256()
    count = 0
    for value in values:
        encoded = json.dumps(
            _normalise_scalar(value), ensure_ascii=False, separators=(",", ":")
        ).encode("utf-8")
        digest.update(len(encoded).to_bytes(8, "little"))
        digest.update(encoded)
        count += 1
    return {"count": count, "sha256": digest.hexdigest()}


def _canonical_array_bytes(values):
    """Return stable bytes for a numeric array under scalar equality semantics."""
    import numpy as np

    array = np.asarray(values)
    if array.dtype.kind not in "biufcSU":
        return None

    # Python/numpy equality treats signed zeroes as equal, and the comparator's
    # scalar contract treats every NaN as the same missing numeric value. Make
    # those cases byte-canonical before hashing.
    if array.dtype.kind in "fc":
        array = array.copy()
        array[array == 0] = 0
        if array.dtype.kind == "f":
            array[np.isnan(array)] = np.nan
        else:
            real = array.real
            imag = array.imag
            real[np.isnan(real)] = np.nan
            imag[np.isnan(imag)] = np.nan

    dtype = array.dtype
    if dtype.byteorder not in ("|", "<"):
        array = array.astype(dtype.newbyteorder("<"), copy=False)
    return np.ascontiguousarray(array).tobytes()


def _update_array_digest(digest, values):
    """Add an array to ``digest`` with an unambiguous length-delimited encoding."""
    import numpy as np

    raw = _canonical_array_bytes(values)
    if raw is None:
        signature = _digest_scalar_sequence(np.asarray(values).flat)
        raw = json.dumps(signature, sort_keys=True, separators=(",", ":")).encode(
            "utf-8"
        )
    digest.update(len(raw).to_bytes(8, "little"))
    digest.update(raw)


def _matrix_chunk_rows(value):
    """Choose a row count whose dense form is approximately four MiB."""
    import numpy as np

    n_rows, n_cols = (int(size) for size in value.shape)
    itemsize = max(1, np.dtype(value.dtype).itemsize)
    bytes_per_row = max(1, n_cols * itemsize)
    return max(1, min(n_rows or 1, _MATRIX_CHUNK_TARGET_BYTES // bytes_per_row))


def _canonical_matrix_block(value, start, stop, require_finite=False):
    """Return a canonical CSR slice, copying at most one bounded row chunk."""
    import numpy as np
    import scipy.sparse as sp

    source = value[start:stop, :]
    source_kind = np.dtype(source.dtype).kind
    if require_finite:
        if source_kind not in "biuf":
            raise ValueError(
                "AnnData X must contain real numeric values, "
                f"found dtype {source.dtype}"
            )
        source_values = source.data if sp.issparse(source) else np.asarray(source)
        if source_kind == "f" and not np.isfinite(source_values).all():
            raise ValueError("AnnData X contains a non-finite value")
        if np.any(source_values < 0):
            raise ValueError("AnnData X contains a negative count value")

    if sp.issparse(source):
        block = source.tocsr(copy=True)
    else:
        block = sp.csr_matrix(np.asarray(source))
    block.sort_indices()
    if require_finite and block.data.dtype.kind in "iu" \
            and not block.has_canonical_format:
        limits = np.iinfo(block.data.dtype)
        for row in range(block.shape[0]):
            first = block.indptr[row]
            last = block.indptr[row + 1]
            indices = block.indices[first:last]
            values = block.data[first:last]
            offset = 0
            while offset < len(indices):
                end = offset + 1
                while end < len(indices) and indices[end] == indices[offset]:
                    end += 1
                if end - offset > 1:
                    total = sum(int(item) for item in values[offset:end])
                    if total < limits.min or total > limits.max:
                        raise ValueError(
                            "integer overflow while canonicalising duplicate "
                            f"AnnData X coordinate ({start + row}, {indices[offset]})"
                        )
                offset = end
    block.sum_duplicates()
    if require_finite and block.data.dtype.kind == "f" \
            and not np.isfinite(block.data).all():
        raise ValueError("AnnData X canonicalisation produced a non-finite value")
    if require_finite and np.any(block.data < 0):
        raise ValueError(
            "AnnData X canonicalisation produced a negative count value"
        )
    block.eliminate_zeros()
    return block


def _matrix_signature(value, require_finite=False):
    """Describe array content using bounded-memory, storage-neutral hashing.

    AnnData count matrices can contain hundreds of millions of non-zero entries.
    Never materialise a COO copy or a Python tuple per entry: backed matrices are
    sliced by row and each small slice is canonicalised as CSR before hashing.
    """
    import numpy as np

    shape = tuple(int(size) for size in value.shape)
    dtype = str(value.dtype)
    digest = hashlib.sha256()
    digest.update(json.dumps({"shape": shape, "dtype": dtype}).encode("utf-8"))

    if len(shape) != 2 or np.dtype(value.dtype).kind not in "biufc":
        if require_finite:
            raise ValueError("AnnData X must be a two-dimensional real numeric matrix")
        array = np.asarray(value)
        itemsize = max(1, array.dtype.itemsize)
        chunk_elements = max(1, _MATRIX_CHUNK_TARGET_BYTES // itemsize)
        flat = array.reshape(-1)
        for start in range(0, flat.size, chunk_elements):
            _update_array_digest(digest, flat[start:start + chunk_elements])
        return {"shape": shape, "dtype": dtype, "sha256": digest.hexdigest()}

    nnz = 0
    chunk_rows = _matrix_chunk_rows(value)
    for start in range(0, shape[0], chunk_rows):
        stop = min(shape[0], start + chunk_rows)
        block = _canonical_matrix_block(
            value, start, stop, require_finite=require_finite
        )
        counts = np.diff(block.indptr)
        rows = np.repeat(np.arange(start, stop, dtype=np.int64), counts)
        columns = block.indices.astype(np.int64, copy=False)
        _update_array_digest(digest, rows)
        _update_array_digest(digest, columns)
        _update_array_digest(digest, block.data)
        nnz += int(block.nnz)

    return {
        "shape": shape,
        "dtype": dtype,
        "nnz": nnz,
        "sha256": digest.hexdigest(),
    }


def _index_signature(index):
    return {
        "name": _normalise_scalar(index.name),
        "dtype": str(index.dtype),
        "values": _digest_scalar_sequence(index),
    }


def _series_dtype_signature(series):
    import pandas as pd

    dtype = series.dtype
    if isinstance(dtype, pd.CategoricalDtype):
        return {
            "dtype": "category",
            "ordered": dtype.ordered,
            "categories": tuple(
                _normalise_scalar(v) for v in dtype.categories.tolist()
            ),
        }
    return {"dtype": str(dtype)}


def _frame_signature(frame):
    """Capture dataframe values, labels and extension-dtype metadata."""
    return {
        "index": _index_signature(frame.index),
        "columns": tuple(_normalise_scalar(c) for c in frame.columns.tolist()),
        "dtypes": tuple(
            _series_dtype_signature(frame.iloc[:, i])
            for i in range(frame.shape[1])
        ),
        "columns_data": tuple(
            _digest_scalar_sequence(frame.iloc[:, i])
            for i in range(frame.shape[1])
        ),
    }


def _nested_signature(value):
    """Canonicalise AnnData's auxiliary mappings (including ``uns``)."""
    import numpy as np
    import pandas as pd
    import scipy.sparse as sp

    if isinstance(value, pd.DataFrame):
        return ("dataframe", _frame_signature(value))
    if isinstance(value, pd.Series):
        return (
            "series",
            _index_signature(value.index),
            _series_dtype_signature(value),
            _digest_scalar_sequence(value),
        )
    if sp.issparse(value) or isinstance(value, np.ndarray):
        return ("array", _matrix_signature(value))
    if isinstance(value, dict):
        return (
            "dict",
            tuple(
                (str(key), _nested_signature(item))
                for key, item in sorted(value.items(), key=lambda pair: str(pair[0]))
            ),
        )
    if isinstance(value, (list, tuple)):
        return (type(value).__name__, tuple(_nested_signature(v) for v in value))
    return ("scalar", _normalise_scalar(value))


def _mapping_signature(mapping):
    return {
        str(key): _nested_signature(mapping[key])
        for key in sorted(mapping.keys(), key=str)
    }


def _h5ad_signature(path):
    """Open an H5AD and build a signature of its customer-visible semantics."""
    adata = _anndata.read_h5ad(path, backed="r")
    try:
        if adata.X is None:
            raise ValueError("AnnData object has no X matrix")

        x_signature = _matrix_signature(adata.X, require_finite=True)
        if "nnz" not in x_signature:
            raise ValueError("AnnData X must be a two-dimensional numeric matrix")
        signature = {
            "X": x_signature,
            "var_names": _index_signature(adata.var_names),
            "obs_names": _index_signature(adata.obs_names),
            "obs": _frame_signature(adata.obs),
            "var": _frame_signature(adata.var),
            "layers": _mapping_signature(adata.layers),
            "obsm": _mapping_signature(adata.obsm),
            "varm": _mapping_signature(adata.varm),
            "obsp": _mapping_signature(adata.obsp),
            "varp": _mapping_signature(adata.varp),
            "uns": _mapping_signature(adata.uns),
        }
        if adata.raw is None:
            signature["raw"] = None
        else:
            signature["raw"] = {
                "X": _matrix_signature(adata.raw.X, require_finite=True),
                "var": _frame_signature(adata.raw.var),
                "varm": _mapping_signature(adata.raw.varm),
            }
        return signature
    finally:
        adata.file.close()


def _compare_h5ad(path_a, path_b, max_diffs, envelope_max_flips=None):
    if os.path.basename(path_a).endswith(".empty.h5ad"):
        sizes = {"a": os.path.getsize(path_a), "b": os.path.getsize(path_b)}
        if sizes == {"a": 0, "b": 0}:
            return EQUAL, {"sentinel": "named empty H5AD", "sizes": sizes}
        return DIFFER, {
            "validation_errors": {
                side: "named .empty.h5ad sentinel must be zero bytes"
                for side, size in sizes.items()
                if size != 0
            },
            "sizes": sizes,
        }

    if not _HAVE_ANNDATA:
        return DIFFER, {
            "dependency_error": "anndata is required to validate .h5ad files"
        }

    signatures = {}
    validation_errors = {}
    for side, path in (("a", path_a), ("b", path_b)):
        try:
            signatures[side] = _h5ad_signature(path)
        except Exception as exc:
            validation_errors[side] = f"{type(exc).__name__}: {exc}"
    if validation_errors:
        return DIFFER, {"validation_errors": validation_errors}

    sig_a = signatures["a"]
    sig_b = signatures["b"]

    mismatches = {}
    x_a = sig_a["X"]
    x_b = sig_b["X"]
    if x_a["shape"] != x_b["shape"]:
        mismatches["shape"] = {"a": x_a["shape"], "b": x_b["shape"]}
    if x_a["dtype"] != x_b["dtype"]:
        mismatches["dtype"] = {"a": x_a["dtype"], "b": x_b["dtype"]}
    if x_a != x_b:
        mismatches["X"] = {
            "a_nnz": x_a["nnz"],
            "b_nnz": x_b["nnz"],
            "reason": "canonical matrix content differs",
        }
    if sig_a["var_names"] != sig_b["var_names"]:
        mismatches["var_names"] = "values or order differ"
    if sig_a["obs_names"] != sig_b["obs_names"]:
        mismatches["obs_names"] = "values or order differ"
    for key in (
        "obs", "var", "layers", "obsm", "varm", "obsp", "varp", "uns", "raw"
    ):
        if sig_a[key] != sig_b[key]:
            mismatches[key] = f"AnnData {key} content differs"

    if not mismatches:
        return EQUAL, {"shape": list(x_a["shape"]), "nnz": x_a["nnz"]}

    if envelope_max_flips is None:
        return DIFFER, {"mismatches": mismatches}

    # Envelope mode. adata.X is obs(barcodes)-by-var(genes) -- the TRANSPOSE of
    # the mtx (genes-by-barcodes). So the per-barcode total == the per-OBS sum
    # (sum over var, axis=1); that is the "column sum" / per-barcode invariant
    # the ambiguity envelope must preserve. The envelope only applies when the
    # difference is confined to X entries: any shape / name-axis mismatch is a
    # structural regression that always DIFFERs.
    structural = {k: v for k, v in mismatches.items() if k != "X"}
    if structural:
        return DIFFER, {"mismatches": mismatches}

    try:
        ok, env_detail = _h5ad_envelope_check(
            path_a, path_b, sig_a, sig_b, envelope_max_flips
        )
    except (OSError, ValueError, OverflowError) as exc:
        return DIFFER, {
            "mismatches": mismatches,
            "validation_errors": {"envelope": f"{type(exc).__name__}: {exc}"},
        }
    detail = {"mismatches": mismatches}
    detail.update(env_detail)
    return (ENVELOPE_OK if ok else DIFFER), detail


def _h5ad_envelope_check(path_a, path_b, sig_a, sig_b, envelope_max_flips):
    """Decide whether the X difference is inside the ambiguity envelope.

    Returns (passes, detail). Passes IFF per-barcode (per-obs) sums are identical
    AND the differing-entry count is within the flip budget.
    """
    import numpy as np

    adata_a = _anndata.read_h5ad(path_a, backed="r")
    adata_b = _anndata.read_h5ad(path_b, backed="r")
    try:
        chunk_rows = min(
            _matrix_chunk_rows(adata_a.X), _matrix_chunk_rows(adata_b.X)
        )
        n_only_a = 0
        n_only_b = 0
        n_diff_cols = 0
        for start in range(0, adata_a.n_obs, chunk_rows):
            stop = min(adata_a.n_obs, start + chunk_rows)
            block_a = _canonical_matrix_block(
                adata_a.X, start, stop, require_finite=True
            )
            block_b = _canonical_matrix_block(
                adata_b.X, start, stop, require_finite=True
            )

            rows_a = np.repeat(
                np.arange(start, stop, dtype=np.int64), np.diff(block_a.indptr)
            )
            rows_b = np.repeat(
                np.arange(start, stop, dtype=np.int64), np.diff(block_b.indptr)
            )
            coordinates_a = rows_a * adata_a.n_vars + block_a.indices
            coordinates_b = rows_b * adata_b.n_vars + block_b.indices
            _shared, indices_a, indices_b = np.intersect1d(
                coordinates_a,
                coordinates_b,
                assume_unique=True,
                return_indices=True,
            )
            values_a = block_a.data[indices_a]
            values_b = block_b.data[indices_b]
            equal_values = values_a == values_b
            n_equal_entries = int(np.count_nonzero(equal_values))
            n_only_a += int(block_a.nnz) - n_equal_entries
            n_only_b += int(block_b.nnz) - n_equal_entries

            for local_row in range(stop - start):
                sums = []
                for side, block in (("a", block_a), ("b", block_b)):
                    values = block.data[
                        block.indptr[local_row]:block.indptr[local_row + 1]
                    ]
                    if values.dtype.kind in "biu":
                        total = sum(int(value) for value in values)
                    elif values.dtype.kind == "f":
                        try:
                            total = math.fsum(float(value) for value in values)
                        except OverflowError as exc:
                            raise ValueError(
                                "numeric overflow while summing AnnData X row "
                                f"{start + local_row} ({side})"
                            ) from exc
                        if not math.isfinite(total):
                            raise ValueError(
                                "non-finite sum for AnnData X row "
                                f"{start + local_row} ({side})"
                            )
                    else:  # guarded by require_finite, retained as fail-closed
                        raise ValueError(
                            f"unsupported AnnData X dtype {values.dtype}"
                        )
                    sums.append(total)
                n_diff_cols += sums[0] != sums[1]
    finally:
        adata_a.file.close()
        adata_b.file.close()

    col_sums_preserved = n_diff_cols == 0
    n_diff_entries = max(n_only_a, n_only_b)

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
    completed = subprocess.run(
        ["samtools"] + args,
        check=False,
        capture_output=True,
        text=True,
    )
    if completed.returncode != 0:
        message = completed.stderr.strip() or completed.stdout.strip()
        raise ValueError(message or f"samtools exited {completed.returncode}")
    return completed.stdout


def _canonical_bam_header(path):
    """Return stable SAM header semantics.

    ``@PG`` records describe the program invocation chain rather than the
    alignments, and ``UR`` is commonly a run-specific reference path. All
    reference identities/lengths and other header records remain part of the
    comparison contract.
    """
    header = []
    sort_orders = []
    for line in _samtools(["view", "-H", path]).splitlines():
        fields = line.split("\t")
        record_type = fields[0]
        if record_type == "@PG":
            continue
        if record_type == "@CO":
            header.append((record_type, "\t".join(fields[1:])))
            continue
        tags = fields[1:]
        if record_type == "@HD":
            sort_order = [tag[3:] for tag in tags if tag.startswith("SO:")]
            if len(sort_order) != 1:
                raise ValueError("SAM @HD must contain exactly one SO tag")
            sort_orders.extend(sort_order)
        if record_type == "@SQ":
            tags = [tag for tag in tags if not tag.startswith("UR:")]
        header.append((record_type, tuple(sorted(tags))))
    if len(sort_orders) != 1:
        raise ValueError("SAM header must contain exactly one @HD record")
    if sort_orders[0] not in {"coordinate", "unsorted", "queryname", "unknown"}:
        raise ValueError(f"unsupported SAM sort order {sort_orders[0]!r}")
    return tuple(header), sort_orders[0]


def _parse_sam_tags(fields):
    tags = {}
    for field in fields:
        parts = field.split(":", 2)
        if len(parts) != 3:
            raise ValueError(f"malformed SAM optional field: {field!r}")
        tag, value_type, value = parts
        if tag in tags:
            raise ValueError(f"duplicate SAM optional tag: {tag}")
        tags[tag] = (value_type, value)
    return tags


def _qname_barcode_umi(qname, tags):
    """Read cell/UMI identity from tags or the pipeline's read-name suffix.

    The current chemistry deduplicates on cell plus genomic SSS start and does
    not encode a separate UMI, so a missing UMI is the canonical empty value.
    """
    suffix = re.search(r"_([ACGTN]+)_([ACGTN]*)$", qname, re.IGNORECASE)
    qname_cell = suffix.group(1) if suffix else None
    qname_umi = suffix.group(2) if suffix else None
    for tag in ("CB", "UB"):
        if tag in tags and tags[tag][0] != "Z":
            raise ValueError(f"SAM {tag} tag must have type Z, found {tags[tag][0]!r}")
    cell = tags.get("CB", (None, qname_cell))[1]
    umi = tags.get("UB", (None, qname_umi))[1]
    if not cell:
        raise ValueError(
            "SAM record has no non-empty cell barcode in CB or the pipeline's "
            "_<cell>_<umi> read-name suffix"
        )
    if umi is None:
        umi = ""
    return cell, umi


def _stable_bam_molecule(line, require_xt=False):
    """Return (coordinate tie, stable molecule identity) for one SAM record.

    Parallel ``umi_tools`` deduplication may select a different source read for
    the same molecule. Its QNAME prefix, sequence, quality, MAPQ, CIGAR, mate
    fields and incidental alignment tags are therefore deliberately excluded.
    The pipeline guarantees the reference/start/strand, cell, UMI and complete
    featureCounts gene assignment (XT) for an equivalent representative.
    """
    fields = line.rstrip("\n").split("\t")
    if len(fields) < 11:
        raise ValueError(f"SAM record has {len(fields)} fields; expected at least 11")

    try:
        flag = int(fields[1])
        position = int(fields[3])
    except ValueError as exc:
        raise ValueError("SAM FLAG and POS must be integers") from exc
    tags = _parse_sam_tags(fields[11:])
    cell, umi = _qname_barcode_umi(fields[0], tags)
    gene = ""
    if "XT" in tags:
        xt_type, gene = tags["XT"]
        if xt_type != "Z":
            raise ValueError(f"SAM XT tag must have type Z, found {xt_type!r}")
        if not gene:
            raise ValueError("SAM XT gene-assignment tag must not be empty")
    elif require_xt:
        raise ValueError(
            "SAM record in a high-confidence annotated/deduplicated BAM "
            "must contain an XT gene-assignment tag"
        )

    coordinate = (fields[2], position)
    molecule = (fields[2], position, bool(flag & 0x10), cell, umi, gene)
    return coordinate, molecule


def _iter_bam_tie_groups(path, require_xt=False):
    """Yield coordinate-tie Counters, bounding memory by the largest tie group."""
    from collections import Counter

    proc = subprocess.Popen(
        ["samtools", "view", path],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    try:
        assert proc.stdout is not None
        coordinate = None
        molecules = Counter()
        for line in proc.stdout:
            next_coordinate, molecule = _stable_bam_molecule(
                line, require_xt=require_xt
            )
            if coordinate is not None and next_coordinate != coordinate:
                yield coordinate, molecules
                molecules = Counter()
            coordinate = next_coordinate
            molecules[molecule] += 1
        if coordinate is not None:
            yield coordinate, molecules
        assert proc.stderr is not None
        stderr = proc.stderr.read().strip()
        returncode = proc.wait()
        if returncode != 0:
            raise ValueError(stderr or f"samtools view exited {returncode}")
    finally:
        if proc.poll() is None:
            proc.kill()
            proc.wait()


def _coordinate_ordered_bam(path, sort_order):
    """Return a coordinate-ordered BAM path and optional temporary-dir owner.

    STAR's configured published BAM is deliberately `SO:unsorted`. Normalise
    such inputs with samtools' external merge sort, capped at 64 MiB per thread,
    rather than retaining every record key in Python memory.
    """
    if sort_order == "coordinate":
        return path, None

    scratch = tempfile.TemporaryDirectory(prefix="compare-bam-")
    output = os.path.join(scratch.name, "coordinate.bam")
    try:
        _samtools([
            "sort", "--no-PG", "-@", "1", "-m", "64M",
            "-o", output, str(path),
        ])
    except Exception:
        scratch.cleanup()
        raise
    return output, scratch


def _bam_requires_xt(path):
    """Whether this configured published BAM contains only assigned records."""
    norm = os.path.normpath(str(path)).replace(os.sep, "/")
    base = os.path.basename(norm)
    return "/deduplication/" in f"/{norm.lstrip('/')}" \
        or base.endswith(".mapped.sorted.filtered.annotated.bam")


def _bam_signature(path, require_xt=False):
    # quickcheck catches missing headers and absent/truncated BGZF EOF blocks;
    # streaming through gzip verifies every BGZF block checksum; consuming
    # every record then validates the decompressed BAM structure and content.
    _samtools(["quickcheck", "-v", path])
    with gzip.open(path, "rb") as bam_stream:
        for _chunk in iter(lambda: bam_stream.read(1 << 20), b""):
            pass
    header, sort_order = _canonical_bam_header(path)
    digest = hashlib.sha256()
    count = 0
    ordered_path, scratch = _coordinate_ordered_bam(path, sort_order)
    try:
        for coordinate, molecules in _iter_bam_tie_groups(
            ordered_path, require_xt=require_xt
        ):
            encoded_coordinate = json.dumps(
                coordinate, ensure_ascii=False, separators=(",", ":")
            ).encode("utf-8")
            digest.update(len(encoded_coordinate).to_bytes(8, "little"))
            digest.update(encoded_coordinate)
            for molecule, multiplicity in sorted(molecules.items()):
                encoded = json.dumps(
                    [molecule, multiplicity],
                    ensure_ascii=False,
                    separators=(",", ":"),
                ).encode("utf-8")
                digest.update(len(encoded).to_bytes(8, "little"))
                digest.update(encoded)
                count += multiplicity
    finally:
        if scratch is not None:
            scratch.cleanup()
    return {
        "header": header,
        "sort_order": sort_order,
        "count": count,
        "sha256": digest.hexdigest(),
    }


def _summarise_counter_difference(counter_a, counter_b, max_diffs):
    only_a = []
    only_b = []
    if max_diffs == 0:
        return only_a, only_b
    for molecule, count in sorted((counter_a - counter_b).items()):
        only_a.append({"molecule": list(molecule), "count": count})
        if len(only_a) >= max_diffs:
            break
    for molecule, count in sorted((counter_b - counter_a).items()):
        only_b.append({"molecule": list(molecule), "count": count})
        if len(only_b) >= max_diffs:
            break
    return only_a, only_b


def _bam_record_diffs(
    path_a,
    path_b,
    max_diffs,
    sort_order_a,
    sort_order_b,
    require_xt_a=False,
    require_xt_b=False,
):
    if max_diffs == 0:
        return []
    missing = object()
    diffs = []
    ordered_a, scratch_a = _coordinate_ordered_bam(path_a, sort_order_a)
    try:
        ordered_b, scratch_b = _coordinate_ordered_bam(path_b, sort_order_b)
    except Exception:
        if scratch_a is not None:
            scratch_a.cleanup()
        raise
    iter_a = _iter_bam_tie_groups(ordered_a, require_xt=require_xt_a)
    iter_b = _iter_bam_tie_groups(ordered_b, require_xt=require_xt_b)
    try:
        groups = itertools.zip_longest(iter_a, iter_b, fillvalue=missing)
        for group_a, group_b in groups:
            if group_a != group_b:
                if group_a is missing or group_b is missing:
                    diffs.append({
                        "coordinate_a": (
                            "<no group>" if group_a is missing else list(group_a[0])
                        ),
                        "coordinate_b": (
                            "<no group>" if group_b is missing else list(group_b[0])
                        ),
                    })
                elif group_a[0] != group_b[0]:
                    diffs.append({
                        "coordinate_a": list(group_a[0]),
                        "coordinate_b": list(group_b[0]),
                    })
                else:
                    only_a, only_b = _summarise_counter_difference(
                        group_a[1], group_b[1], max_diffs
                    )
                    diffs.append({
                        "coordinate": list(group_a[0]),
                        "only_a": only_a,
                        "only_b": only_b,
                    })
                if len(diffs) >= max_diffs:
                    break
    finally:
        iter_a.close()
        iter_b.close()
        if scratch_a is not None:
            scratch_a.cleanup()
        if scratch_b is not None:
            scratch_b.cleanup()
    return diffs


def _compare_bam(path_a, path_b, max_diffs):
    if not _HAVE_SAMTOOLS:
        return DIFFER, {
            "dependency_error": "samtools is required to validate .bam files"
        }

    signatures = {}
    validation_errors = {}
    require_xt = {
        "a": _bam_requires_xt(path_a),
        "b": _bam_requires_xt(path_b),
    }
    for side, path in (("a", path_a), ("b", path_b)):
        try:
            signatures[side] = _bam_signature(
                path, require_xt=require_xt[side]
            )
        except Exception as exc:
            validation_errors[side] = f"{type(exc).__name__}: {exc}"
    if validation_errors:
        return DIFFER, {"validation_errors": validation_errors}

    sig_a = signatures["a"]
    sig_b = signatures["b"]
    if sig_a == sig_b:
        return EQUAL, {"count": sig_a["count"]}

    detail = {}
    if sig_a["header"] != sig_b["header"]:
        detail["header"] = "stable SAM header records differ"
    if sig_a["count"] != sig_b["count"]:
        detail["count"] = {"a": sig_a["count"], "b": sig_b["count"]}
    if sig_a["sha256"] != sig_b["sha256"]:
        detail["record_sha256"] = {
            "a": sig_a["sha256"],
            "b": sig_b["sha256"],
        }
        detail["record_diffs"] = _bam_record_diffs(
            path_a,
            path_b,
            max_diffs,
            sig_a["sort_order"],
            sig_b["sort_order"],
            require_xt["a"],
            require_xt["b"],
        )
    return DIFFER, detail


# ----------------------------------------------------------------------------
# HTML
# ----------------------------------------------------------------------------

_HTML_VOID_ELEMENTS = {
    "area", "base", "br", "col", "embed", "hr", "img", "input", "link",
    "meta", "param", "source", "track", "wbr",
}

_HTML_OPTIONAL_END_ELEMENTS = {
    "html", "head", "body", "p", "li", "dt", "dd", "rt", "rp",
    "optgroup", "option", "colgroup", "thead", "tbody", "tfoot", "tr",
    "td", "th",
}

_HTML_AUTO_CLOSE_ON_SAME = {
    "li", "dt", "dd", "rt", "rp", "optgroup", "option", "thead", "tbody",
    "tfoot", "tr", "td", "th",
}

_HTML_P_CLOSING_START_TAGS = {
    "address", "article", "aside", "blockquote", "details", "div", "dl",
    "fieldset", "figcaption", "figure", "footer", "form", "h1", "h2", "h3",
    "h4", "h5", "h6", "header", "hgroup", "hr", "main", "menu", "nav",
    "ol", "p", "pre", "search", "section", "table", "ul",
}

_HTML_INTRINSIC_CONTENT_ELEMENTS = {
    "audio", "canvas", "embed", "iframe", "img", "input", "object", "svg",
    "video",
}

_EMPTY_HTML_SENTINEL_SUFFIXES = (
    "_counts_pdf_with_threshold.html",
    "_barnyard_plot.html",
    "_pdf_with_cutoff.html",
)


class _HTMLValidator(html.parser.HTMLParser):
    """Validate structure and retain parsed elements needed by output contracts."""

    def __init__(self):
        super().__init__(convert_charrefs=False)
        self.stack = []
        self.error = None
        self.start_tags = 0
        self.has_content = False
        self.elements = []
        self.text_by_tag = {}
        self.scripts = []
        self._captures = []

    def _set_error(self, message):
        if self.error is None:
            line, column = self.getpos()
            self.error = f"line {line}, column {column}: {message}"

    def handle_starttag(self, tag, attrs):
        tag = tag.lower()
        self.start_tags += 1
        attributes = {
            str(name).lower(): "" if value is None else value
            for name, value in attrs
        }
        self.elements.append({
            "tag": tag,
            "attrs": attributes,
            "ancestors": tuple(self.stack),
        })
        if tag in _HTML_INTRINSIC_CONTENT_ELEMENTS:
            self.has_content = True
        if self.stack and self.stack[-1] == "p" and tag in _HTML_P_CLOSING_START_TAGS:
            self.stack.pop()
        if self.stack and self.stack[-1] == tag and tag in _HTML_AUTO_CLOSE_ON_SAME:
            self.stack.pop()
        if tag not in _HTML_VOID_ELEMENTS:
            self.stack.append(tag)
        if tag in {"title", "h1", "h2", "h3", "pre", "script"}:
            self._captures.append({"tag": tag, "attrs": attributes, "parts": []})

    def handle_startendtag(self, tag, attrs):
        tag = tag.lower()
        self.start_tags += 1
        attributes = {
            str(name).lower(): "" if value is None else value
            for name, value in attrs
        }
        self.elements.append({
            "tag": tag,
            "attrs": attributes,
            "ancestors": tuple(self.stack),
        })
        if tag in _HTML_INTRINSIC_CONTENT_ELEMENTS:
            self.has_content = True

    def handle_endtag(self, tag):
        tag = tag.lower()
        if tag in _HTML_VOID_ELEMENTS:
            self._set_error(f"void element </{tag}> must not have an end tag")
            return
        if tag not in self.stack:
            self._set_error(f"unexpected closing tag </{tag}>")
            return
        index = len(self.stack) - 1 - self.stack[::-1].index(tag)
        unclosed = self.stack[index + 1:]
        if any(open_tag not in _HTML_OPTIONAL_END_ELEMENTS for open_tag in unclosed):
            self._set_error(
                f"closing </{tag}> before <{self.stack[-1]}> was closed"
            )
            return
        for capture_index in range(len(self._captures) - 1, -1, -1):
            capture = self._captures[capture_index]
            if capture["tag"] == tag:
                text = "".join(capture["parts"])
                if tag == "script":
                    self.scripts.append({"attrs": capture["attrs"], "body": text})
                else:
                    self.text_by_tag.setdefault(tag, []).append(text)
                del self._captures[capture_index]
                break
        del self.stack[index:]

    def handle_data(self, data):
        if data.strip():
            self.has_content = True
        for capture in self._captures:
            capture["parts"].append(data)

    def handle_entityref(self, name):
        self.has_content = True
        for capture in self._captures:
            capture["parts"].append(f"&{name};")

    def handle_charref(self, name):
        self.has_content = True
        for capture in self._captures:
            capture["parts"].append(f"&#{name};")

    def finish(self):
        if self.error:
            return self.error
        required_end_tags = [
            tag for tag in self.stack if tag not in _HTML_OPTIONAL_END_ELEMENTS
        ]
        if required_end_tags:
            return f"unclosed tag <{required_end_tags[-1]}>"
        if self.start_tags == 0:
            return "no HTML elements"
        if not self.has_content:
            return "HTML contains no renderable content"
        return None


def _parse_html(path):
    try:
        raw = _read_bytes(path)
        text = raw.decode("utf-8", "strict")
    except UnicodeDecodeError as exc:
        return None, f"invalid UTF-8: {exc}"
    if not text.strip():
        return None, "empty HTML"

    validator = _HTMLValidator()
    try:
        validator.feed(text)
        validator.close()
    except Exception as exc:
        return None, f"HTML parse error: {type(exc).__name__}: {exc}"
    error = validator.finish()
    return (None, error) if error else (validator, None)


def _normalised_html_text(document, tag):
    return [
        " ".join(html.unescape(value).split())
        for value in document.text_by_tag.get(tag, [])
    ]


def _has_html_element(document, tag=None, element_id=None, class_token=None):
    for element in document.elements:
        attrs = element["attrs"]
        if tag is not None and element["tag"] != tag:
            continue
        if element_id is not None and attrs.get("id") != element_id:
            continue
        if class_token is not None and class_token not in attrs.get("class", "").split():
            continue
        return True
    return False


def _validate_multiqc_html(document):
    titles = _normalised_html_text(document, "title")
    if not any("multiqc" in title.lower() for title in titles):
        return "missing parsed MultiQC title"
    config_scripts = [
        script for script in document.scripts
        if script["attrs"].get("id") == "mqc_config"
    ]
    if not config_scripts:
        return "missing parsed MultiQC configuration script"
    if len(config_scripts) != 1:
        return "MultiQC report must contain exactly one configuration script"
    config_script = config_scripts[0]
    config_type = config_script["attrs"].get("type", "").strip().lower()
    if config_type != "application/json":
        return "MultiQC configuration script must have type application/json"
    try:
        config = _load_json_strict(config_script["body"])
    except ValueError as exc:
        return f"invalid MultiQC configuration JSON: {exc}"
    if not isinstance(config, dict) or not config:
        return "MultiQC configuration JSON must be a non-empty object"
    required_config = {
        "decimalPoint_format",
        "num_datasets_plot_limit",
        "sample_names_rename",
        "show_hide_mode",
        "show_hide_patterns",
        "show_hide_regex",
        "thousandsSep_format",
    }
    missing_config = sorted(required_config - config.keys())
    if missing_config:
        return f"MultiQC configuration lacks pinned 1.14 fields: {missing_config}"
    for key in ("decimalPoint_format", "thousandsSep_format"):
        if config[key] is not None and (
            not isinstance(config[key], str) or not config[key]
        ):
            return f"MultiQC configuration {key} has an invalid value"
    plot_limit = config["num_datasets_plot_limit"]
    if isinstance(plot_limit, bool) or not isinstance(plot_limit, int) \
            or plot_limit <= 0:
        return "MultiQC configuration num_datasets_plot_limit must be positive"
    for key in (
        "sample_names_rename", "show_hide_mode", "show_hide_patterns",
        "show_hide_regex",
    ):
        if not isinstance(config[key], list):
            return f"MultiQC configuration {key} must be a list"
    if not _has_html_element(document, "h1", "page_title"):
        return "missing parsed MultiQC report heading"
    has_module = any(
        element["tag"] in {"div", "section"}
        and element["attrs"].get("id", "").startswith("mqc-module-section-")
        for element in document.elements
    )
    if not has_module:
        return "missing parsed MultiQC module section"
    return None


def _is_executable_inline_script(script):
    """Whether an HTML script body executes under the normal MIME rules."""
    attrs = script["attrs"]
    if "src" in attrs:
        return False
    mime = attrs.get("type", "").strip().lower().split(";", 1)[0]
    return mime in {
        "",
        "module",
        "text/javascript",
        "application/javascript",
        "text/ecmascript",
        "application/ecmascript",
    }


def _javascript_control_before_closing_parenthesis(masked, closing):
    """Return the control keyword whose condition ends at ``closing``, if any."""
    depth = 0
    cursor = closing
    while cursor >= 0:
        char = masked[cursor]
        if char == ")":
            depth += 1
        elif char == "(":
            depth -= 1
            if depth == 0:
                cursor -= 1
                while cursor >= 0 and masked[cursor].isspace():
                    cursor -= 1
                end = cursor + 1
                while cursor >= 0 and (
                    masked[cursor].isalnum() or masked[cursor] in "_$"
                ):
                    cursor -= 1
                return "".join(masked[cursor + 1:end])
        cursor -= 1
    return ""


def _javascript_code_mask(source):
    """Return JavaScript with literal/comment contents blanked in place.

    This is a deliberately small lexical pass, not a JavaScript evaluator. It
    preserves positions and string delimiters while replacing the contents of
    quoted/template literals and comments with spaces. Semantic landmark
    searches can therefore distinguish executable tokens from marker text in a
    string or comment without executing customer-produced HTML.
    """
    masked = list(source)
    length = len(source)
    index = 0
    while index < length:
        char = source[index]
        if char in {"'", '"', "`"}:
            quote = char
            index += 1
            while index < length:
                char = source[index]
                if char == "\\":
                    masked[index] = " "
                    index += 1
                    if index < length:
                        masked[index] = " "
                        index += 1
                    continue
                if char == quote:
                    index += 1
                    break
                if quote != "`" and char in "\r\n":
                    raise ValueError("unterminated JavaScript string literal")
                masked[index] = "\n" if char == "\n" else " "
                index += 1
            else:
                raise ValueError("unterminated JavaScript literal")
            continue
        if source.startswith("//", index) or source.startswith("<!--", index):
            while index < length and source[index] not in "\r\n":
                masked[index] = " "
                index += 1
            continue
        if source.startswith("/*", index):
            end = source.find("*/", index + 2)
            if end < 0:
                raise ValueError("unterminated JavaScript block comment")
            for offset in range(index, end + 2):
                if source[offset] not in "\r\n":
                    masked[offset] = " "
            index = end + 2
            continue
        if char == "/":
            # JavaScript regex literals can contain quote/comment characters.
            # Distinguish the common expression-prefix contexts from division,
            # then blank a complete literal so those characters cannot confuse
            # the string/comment pass or act as semantic marker decoys.
            previous = index - 1
            while previous >= 0 and masked[previous].isspace():
                previous -= 1
            can_start_regex = previous < 0 or masked[previous] in "([{:;,=!?&|+-*%^~<>"
            if previous >= 0 and masked[previous] == ")":
                can_start_regex = (
                    _javascript_control_before_closing_parenthesis(masked, previous)
                    in {"catch", "for", "if", "switch", "while", "with"}
                )
            if not can_start_regex and previous >= 0 \
                    and (masked[previous].isalnum() or masked[previous] in "_$"):
                start = previous
                while start >= 0 and (
                    masked[start].isalnum() or masked[start] in "_$"
                ):
                    start -= 1
                can_start_regex = "".join(masked[start + 1:previous + 1]) in {
                    "await", "case", "delete", "do", "else", "in", "instanceof",
                    "new", "of", "return", "throw", "typeof", "void", "yield",
                }
            if can_start_regex:
                cursor = index + 1
                in_character_class = False
                closing = None
                while cursor < length and source[cursor] not in "\r\n":
                    if source[cursor] == "\\":
                        cursor += 2
                        continue
                    if source[cursor] == "[":
                        in_character_class = True
                    elif source[cursor] == "]":
                        in_character_class = False
                    elif source[cursor] == "/" and not in_character_class:
                        closing = cursor
                        break
                    cursor += 1
                if closing is not None:
                    for offset in range(index + 1, closing):
                        masked[offset] = " "
                    index = closing + 1
                    while index < length and source[index].isalpha():
                        index += 1
                    continue
        index += 1
    return "".join(masked)


def _position_between_slash_delimiters(code, start, end):
    """Fail closed when a marker lies inside a same-line regex-like lexeme.

    The main lexer identifies regex literals from expression context so it can
    safely traverse Nextflow's bundled JavaScript. JavaScript's division/regex
    lexical goal is grammar-dependent, however. This target-local guard catches
    any landmark left between slash delimiters even in rare statement contexts
    that a context heuristic cannot infer. Generated pipeline landmark calls
    are never arithmetic expressions between slash operators.
    """
    line_start = code.rfind("\n", 0, start) + 1
    line_end = code.find("\n", end)
    if line_end < 0:
        line_end = len(code)

    slashes = []
    for index in range(line_start, line_end):
        if code[index] != "/":
            continue
        escapes = 0
        previous = index - 1
        while previous >= line_start and code[previous] == "\\":
            escapes += 1
            previous -= 1
        if escapes % 2 == 0:
            slashes.append(index)
    return any(left < start and right >= end
               for left in slashes for right in slashes if left < right)


def _javascript_string_at(source, index):
    """Decode a simple JavaScript string literal beginning at ``index``."""
    if index >= len(source) or source[index] not in {"'", '"'}:
        raise ValueError("expected a quoted JavaScript string argument")
    quote = source[index]
    index += 1
    value = []
    simple_escapes = {
        "b": "\b", "f": "\f", "n": "\n", "r": "\r", "t": "\t",
        "v": "\v", "0": "\0", "\\": "\\", "'": "'", '"': '"',
        "/": "/",
    }
    while index < len(source):
        char = source[index]
        if char == quote:
            return "".join(value), index + 1
        if char in "\r\n":
            raise ValueError("newline in JavaScript string argument")
        if char != "\\":
            value.append(char)
            index += 1
            continue
        index += 1
        if index >= len(source):
            break
        escaped = source[index]
        if escaped in simple_escapes:
            value.append(simple_escapes[escaped])
            index += 1
            continue
        if escaped in {"u", "x"}:
            width = 4 if escaped == "u" else 2
            digits = source[index + 1:index + 1 + width]
            if len(digits) != width or not re.fullmatch(r"[0-9a-fA-F]+", digits):
                raise ValueError("invalid JavaScript hexadecimal string escape")
            value.append(chr(int(digits, 16)))
            index += width + 1
            continue
        raise ValueError(f"unsupported JavaScript string escape \\{escaped}")
    raise ValueError("unterminated JavaScript string argument")


def _plotly_call_targets(document):
    """Return DOM ids passed as the first argument to real Plotly.newPlot calls."""
    targets = set()
    for script in document.scripts:
        if not _is_executable_inline_script(script):
            continue
        source = script["body"]
        if "Plotly" not in source or "newPlot" not in source:
            continue
        code = _javascript_code_mask(source)
        for match in re.finditer(r"\bPlotly\s*\.\s*newPlot\s*\(", code):
            if _position_between_slash_delimiters(
                code, match.start(), match.end()
            ):
                continue
            index = match.end()
            while index < len(code) and code[index].isspace():
                index += 1
            target, _end = _javascript_string_at(source, index)
            if not target:
                raise ValueError("Plotly.newPlot target id must not be empty")
            targets.add(target)
    return targets


_NEXTFLOW_WINDOW_DATA_TRAILERS = {
    "execution_report.html": ";",
    "execution_timeline.html": (
        ";\n\n"
        "  var prot = ((\"https:\" == document.location.protocol) "
        "? \"https://\" : \"http://\");\n"
        "  document.write(unescape(\"%3Clink href='\" + prot + "
        "\"fonts.googleapis.com/css?family=Lato' rel='stylesheet' "
        "type='text/css' %3E%3C/link%3E\"));"
    ),
}


def _nextflow_window_data(document, output_name):
    """Parse Nextflow 26.04.1's ``window.data`` JavaScript object literal.

    The report serializer uses the JavaScript-only ``\\'`` escape inside
    double-quoted command strings. Removing that redundant escape makes the
    object strict JSON. After it, accept only the filename-specific trailer
    emitted by pinned Nextflow; unparsed executable suffixes fail closed.
    """
    decoder = json.JSONDecoder(
        parse_constant=_reject_json_constant,
        object_pairs_hook=_json_object_without_duplicates,
    )
    if output_name not in _NEXTFLOW_WINDOW_DATA_TRAILERS:
        raise ValueError(
            f"no pinned Nextflow window.data trailer for {output_name!r}"
        )

    parsed = []
    candidate_errors = []
    for script in document.scripts:
        if not _is_executable_inline_script(script):
            continue
        if not re.search(r"\bwindow\s*\.\s*data\b", script["body"]):
            continue
        try:
            code = _javascript_code_mask(script["body"])
        except ValueError as exc:
            raise ValueError(f"invalid Nextflow script syntax: {exc}") from exc
        for match in re.finditer(r"\bwindow\s*\.\s*data\s*=\s*", code):
            if _position_between_slash_delimiters(
                code, match.start(), match.end()
            ):
                continue
            candidate = script["body"][match.end():].replace("\\'", "'")
            try:
                value, end = decoder.raw_decode(candidate)
                _validate_finite_json(value)
                if not isinstance(value, dict):
                    raise ValueError("Nextflow window.data must be an object")
                trailer = candidate[end:].strip(" \t\r\n")
                expected_trailer = _NEXTFLOW_WINDOW_DATA_TRAILERS[output_name]
                if trailer != expected_trailer:
                    raise ValueError(
                        f"unexpected executable trailer after Nextflow window.data "
                        f"in {output_name}"
                    )
            except (json.JSONDecodeError, ValueError) as exc:
                candidate_errors.append(str(exc))
                continue
            parsed.append(value)
    if candidate_errors:
        raise ValueError(
            f"invalid Nextflow window.data object: {candidate_errors[0]}"
        )
    if len(parsed) > 1:
        raise ValueError("multiple executable Nextflow window.data assignments")
    if parsed:
        return parsed[0]
    raise ValueError("missing parsed Nextflow window.data script")


_ERROR_PAGE_HEADING = re.compile(
    r"(?:^|\b)(?:404\s+not\s+found|500\s+internal\s+server\s+error|"
    r"internal\s+server\s+error|bad\s+gateway|service\s+unavailable|"
    r"uncaught\s+exception|traceback|failed|failure|error)"
    r"(?:\b|$)",
    re.IGNORECASE,
)


def _validate_execution_report(document):
    titles = _normalised_html_text(document, "title")
    if not any("Nextflow Workflow Report" in title for title in titles):
        return "missing Nextflow workflow-report title"
    if not _has_html_element(document, "nav", "nf-report-navbar"):
        return "missing Nextflow report navigation element"
    if not _has_html_element(document, "table", "tasks_table"):
        return "missing Nextflow task table"
    try:
        data = _nextflow_window_data(document, "execution_report.html")
    except ValueError as exc:
        return str(exc)
    trace = data.get("trace")
    summary = data.get("summary")
    if not isinstance(trace, list) or not trace:
        return "Nextflow report trace must contain task records"
    if not isinstance(summary, list) or not summary:
        return "Nextflow report summary must contain process records"
    required = {"process", "name", "tag", "status", "exit"}
    for index, record in enumerate(trace):
        if not isinstance(record, dict) or not required <= record.keys():
            return f"Nextflow report trace record {index} lacks stable task fields"
        if not all(isinstance(record[field], str) and record[field].strip()
                   for field in ("process", "name")):
            return f"Nextflow report trace record {index} has an empty task identity"
        if record["tag"] is not None and not isinstance(record["tag"], str):
            return f"Nextflow report trace record {index} has an invalid tag"
        exit_value = record["exit"]
        if type(exit_value) is int:
            exit_code = exit_value
        elif isinstance(exit_value, str) and re.fullmatch(r"\d+", exit_value):
            exit_code = int(exit_value)
        else:
            return f"Nextflow report trace record {index} has an invalid exit code"
        status = record["status"]
        if not isinstance(status, str) or status not in _TRACE_SUCCESS_STATUSES \
                or exit_code != 0:
            return f"Nextflow report records unsuccessful task {record['name']!r}"
    if not all(isinstance(record, dict)
               and isinstance(record.get("process"), str)
               and record["process"].strip()
               for record in summary):
        return "Nextflow report summary contains an invalid process record"
    return None


def _validate_execution_timeline(document):
    headings = _normalised_html_text(document, "h3")
    if "Processes execution timeline" not in headings:
        return "missing Nextflow execution-timeline heading"
    if not _has_html_element(document, "div", "timeline"):
        return "missing Nextflow timeline element"
    for element_id in ("label_launch", "label_elapsed", "label_legend"):
        if not _has_html_element(document, "span", element_id):
            return f"missing Nextflow timeline label #{element_id}"
    try:
        data = _nextflow_window_data(document, "execution_timeline.html")
    except ValueError as exc:
        return str(exc)
    if not isinstance(data.get("elapsed"), str) or not data["elapsed"].strip():
        return "Nextflow timeline has no elapsed duration"
    beginning = data.get("beginningMillis")
    ending = data.get("endingMillis")
    if type(beginning) is not int or type(ending) is not int or ending < beginning:
        return "Nextflow timeline has invalid execution bounds"
    processes = data.get("processes")
    if not isinstance(processes, list) or not processes:
        return "Nextflow timeline must contain process records"
    for index, process in enumerate(processes):
        if not isinstance(process, dict) \
                or not isinstance(process.get("label"), str) \
                or not process["label"].strip():
            return f"Nextflow timeline process {index} has no label"
        if not isinstance(process.get("cached"), bool):
            return f"Nextflow timeline process {index} has no cached-state boolean"
        times = process.get("times")
        if not isinstance(times, list) or not times:
            return f"Nextflow timeline process {index} has no timing segments"
        for timing in times:
            start = timing.get("starting_time") if isinstance(timing, dict) else None
            end = timing.get("ending_time") if isinstance(timing, dict) else None
            if type(start) is not int or type(end) is not int or end < start:
                return f"Nextflow timeline process {index} has invalid timing data"
    return None


def _mermaid_diagrams(document):
    pre_elements = [e for e in document.elements if e["tag"] == "pre"]
    return [
        text for element, text in zip(
            pre_elements, document.text_by_tag.get("pre", [])
        )
        if "mermaid" in element["attrs"].get("class", "").split()
    ]


def _validate_pipeline_dag(document):
    diagrams = _mermaid_diagrams(document)
    if not diagrams:
        return "missing parsed Mermaid DAG element"
    diagram = html.unescape(diagrams[0])
    if not re.search(r"(?m)^\s*flowchart\s+(?:TB|TD|BT|RL|LR)\b", diagram):
        return "Mermaid DAG does not declare a flowchart"
    if not re.search(r"(?m)^\s*[A-Za-z_]\w*\s*(?:\[|\()", diagram):
        return "Mermaid DAG contains no workflow nodes"
    has_mermaid_script = False
    for script in document.scripts:
        if script["attrs"].get("type", "").lower() != "module":
            continue
        if not _is_executable_inline_script(script):
            continue
        if "mermaid" not in script["body"] or "initialize" not in script["body"]:
            continue
        try:
            code = _javascript_code_mask(script["body"])
        except ValueError as exc:
            return f"invalid Mermaid module script: {exc}"
        imports = [
            match for match in re.finditer(r"\bimport\s+mermaid\b", code)
            if not _position_between_slash_delimiters(
                code, match.start(), match.end()
            )
        ]
        initialisers = [
            match
            for match in re.finditer(
                r"\bmermaid\s*\.\s*initialize\s*\(", code
            )
            if not _position_between_slash_delimiters(
                code, match.start(), match.end()
            )
        ]
        if imports and initialisers:
            has_mermaid_script = True
            break
    if not has_mermaid_script:
        return "missing parsed Mermaid module initialisation"
    return None


_NEXTFLOW_HTML_NAMES = {
    "execution_report.html",
    "execution_timeline.html",
    "pipeline_dag.html",
}


def _nextflow_html_signature(base, document):
    """Return stable Nextflow semantics while excluding timings/resources/IDs."""
    from collections import Counter

    if base == "execution_report.html":
        data = _nextflow_window_data(document, base)
        records = Counter(
            (
                str(record["process"]),
                str(record["name"]),
                "" if record["tag"] is None else str(record["tag"]),
                str(record["status"]),
                int(record["exit"]),
            )
            for record in data["trace"]
        )
        summaries = Counter(str(record["process"]) for record in data["summary"])
        return (
            "report",
            ("tasks", tuple(sorted(records.items()))),
            ("process_summaries", tuple(sorted(summaries.items()))),
        )

    if base == "execution_timeline.html":
        data = _nextflow_window_data(document, base)
        processes = Counter(
            (str(process["label"]), process["cached"])
            for process in data["processes"]
        )
        return ("processes", tuple(sorted(processes.items())))

    diagram = html.unescape(_mermaid_diagrams(document)[0])
    lines = diagram.splitlines()
    first_flowchart = next(
        index for index, line in enumerate(lines)
        if re.match(r"^\s*flowchart\s+", line)
    )
    # The preceding Mermaid theme block is presentation. The flowchart lines
    # themselves are the workflow topology and labels and must remain stable.
    graph = tuple(line.strip() for line in lines[first_flowchart:] if line.strip())
    return ("flowchart", graph)


def _validate_html_semantics(path, document):
    """Require parsed landmarks belonging to the published output filename."""
    base = os.path.basename(path)

    titles = _normalised_html_text(document, "title")
    h1 = _normalised_html_text(document, "h1")
    h2 = _normalised_html_text(document, "h2")
    # Identify an error *page* from its title/top heading. Report sections may
    # legitimately discuss errors and must not be mistaken for an error page.
    headings = titles + h1[:1] + ([] if h1 else h2[:1])
    if any(_ERROR_PAGE_HEADING.search(heading) for heading in headings):
        return "document identifies itself as an error or failure page"

    if base == "consolidated_report.html":
        if not any(title.startswith("CS Genetics scRNA-seq report") for title in titles):
            return "missing consolidated-report title"
        required_elements = (
            ("main", "main-content", None),
            ("section", "overview", None),
            ("section", "cross-sample", None),
            ("table", None, "cs-xsample-table"),
        )
        missing = [
            element_id or f"{tag}.{class_token}"
            for tag, element_id, class_token in required_elements
            if not _has_html_element(document, tag, element_id, class_token)
        ]
        return None if not missing else f"missing consolidated-report elements: {missing}"

    if base == "multisample_multiqc.html" or base.endswith("_multiqc.html"):
        return _validate_multiqc_html(document)

    if base == "execution_report.html":
        return _validate_execution_report(document)

    if base == "execution_timeline.html":
        return _validate_execution_timeline(document)

    if base == "pipeline_dag.html":
        return _validate_pipeline_dag(document)

    if (
        base == "multisample_qc_cascade.html"
        or base.endswith(".qc_cascade.html")
        or base.endswith("_counts_pdf_with_threshold.html")
        or base.endswith("_barnyard_plot.html")
        or base.endswith("_pdf_with_cutoff.html")
    ):
        graph_ids = {
            element["attrs"].get("id")
            for element in document.elements
            if element["tag"] == "div"
            and "plotly-graph-div" in element["attrs"].get("class", "").split()
            and element["attrs"].get("id")
        }
        try:
            plot_targets = _plotly_call_targets(document)
        except ValueError as exc:
            return f"invalid Plotly script: {exc}"
        if graph_ids & plot_targets:
            return None
        return "missing Plotly.newPlot call targeting a parsed Plotly graph element"

    return "unrecognised published HTML filename; no semantic contract is defined"


def _compare_html(path_a, path_b, max_diffs):
    size_a = os.path.getsize(path_a)
    size_b = os.path.getsize(path_b)
    base = os.path.basename(path_a)
    if base.endswith(_EMPTY_HTML_SENTINEL_SUFFIXES) and size_a == size_b == 0:
        return PRESENT, {
            "size_a": 0,
            "size_b": 0,
            "reason": "both files are an explicitly named empty plot sentinel",
        }

    documents = {}
    validation_errors = {}
    for side, path in (("a", path_a), ("b", path_b)):
        document, error = _parse_html(path)
        if error is None:
            error = _validate_html_semantics(path, document)
        if error:
            validation_errors[side] = error
        else:
            documents[side] = document
    if validation_errors:
        return DIFFER, {
            "size_a": size_a,
            "size_b": size_b,
            "validation_errors": validation_errors,
        }

    if base in _NEXTFLOW_HTML_NAMES:
        signature_a = _nextflow_html_signature(base, documents["a"])
        signature_b = _nextflow_html_signature(base, documents["b"])
        if signature_a == signature_b:
            return EQUAL, {
                "reason": "stable parsed Nextflow semantics are equal",
            }
        encoded_a = json.dumps(signature_a, separators=(",", ":")).encode("utf-8")
        encoded_b = json.dumps(signature_b, separators=(",", ":")).encode("utf-8")
        return DIFFER, {
            "reason": "stable parsed Nextflow semantics differ",
            "semantic_sha256_a": hashlib.sha256(encoded_a).hexdigest(),
            "semantic_sha256_b": hashlib.sha256(encoded_b).hexdigest(),
        }

    # Generated presentation artefacts contain deliberately volatile Plotly
    # element IDs and run provenance. Their structure must be valid, but those
    # values are not numerical pipeline outputs.
    return PRESENT, {
        "size_a": size_a,
        "size_b": size_b,
        "reason": "both files satisfy their structural and output-specific contract",
    }


# ----------------------------------------------------------------------------
# Dispatch
# ----------------------------------------------------------------------------

_COMPARATORS = {
    "TEXT_EXACT": _compare_text_exact,
    "DEDUP_LOG": _compare_dedup,
    "CONFIG": _compare_config,
    "TRACE": _compare_trace,
    "MULTIQC_LOG": _compare_multiqc_log,
    "MULTIQC_JSON": _compare_multiqc_json,
    "MULTIQC_SOURCES": _compare_multiqc_sources,
    "GZ_TEXT_EXACT": _compare_gz_text_exact,
    "MTX": _compare_mtx,
    "H5AD": _compare_h5ad,
    "BAM": _compare_bam,
    "HTML": _compare_html,
    "BINARY_EXACT": _compare_binary_exact,
}


# Classes whose comparator accepts the envelope flip budget. Every other class
# keeps its documented strict contract; the envelope must NOT relax them.
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

def run(
    dir_a,
    dir_b,
    max_diffs,
    envelope_max_flips=None,
    require_complete=True,
):
    files_a = list_relfiles(dir_a)
    files_b = list_relfiles(dir_b)

    only_a = sorted(files_a - files_b)
    only_b = sorted(files_b - files_a)
    missing_both = sorted(
        _REQUIRED_PUBLISHED_FILES - (files_a | files_b)
    ) if require_complete else []
    empty_subset = not require_complete and not files_a and not files_b
    paired = sorted(files_a & files_b)

    results = []
    file_set_mismatch = bool(only_a or only_b or missing_both or empty_subset)

    for rel in only_a:
        results.append({
            "path": rel, "class": "FILE_SET_MISMATCH",
            "verdict": MISSING,
            "detail": {
                "present_in": "a",
                "required": rel in _REQUIRED_PUBLISHED_FILES,
            },
        })
    for rel in only_b:
        results.append({
            "path": rel, "class": "FILE_SET_MISMATCH",
            "verdict": MISSING,
            "detail": {
                "present_in": "b",
                "required": rel in _REQUIRED_PUBLISHED_FILES,
            },
        })
    for rel in missing_both:
        results.append({
            "path": rel, "class": "FILE_SET_MISMATCH",
            "verdict": MISSING,
            "detail": {"missing_from": ["a", "b"], "required": True},
        })
    if empty_subset:
        results.append({
            "path": ".", "class": "FILE_SET_MISMATCH",
            "verdict": MISSING,
            "detail": {
                "missing_from": ["a", "b"],
                "required": "at least one paired output in subset mode",
            },
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
    # verdict (EQUAL/ENVELOPE_OK/PRESENT) never fails; ENVELOPE_OK is the
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
        "require_complete": require_complete,
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
    if r["verdict"] == DIFFER and "validation_errors" in detail:
        lines.append(f"      validation failed: {json.dumps(detail['validation_errors'])}")
        return lines
    if cls == "MTX" and r["verdict"] == DIFFER:
        if "metadata_a" in detail:
            lines.append(
                f"      metadata differs: a={detail['metadata_a']} "
                f"b={detail['metadata_b']}"
            )
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
    elif cls in (
        "H5AD", "CONFIG", "TRACE", "MULTIQC_JSON", "MULTIQC_SOURCES",
        "BAM", "HTML",
    ) \
            and r["verdict"] == DIFFER:
        lines.append(f"      {json.dumps(detail)}")
    elif r["verdict"] == PRESENT and detail.get("warn"):
        lines.append(f"      {detail['warn']}")
    return lines


def print_report(results, summary, max_diffs):
    print(f"Comparing:\n  A = {summary['dir_a']}\n  B = {summary['dir_b']}")
    emf = summary.get("envelope_max_flips")
    mode = "strict" if emf is None else f"envelope (max-flips={emf})"
    scope = "published-tree" if summary.get("require_complete", True) \
        else "explicit-subset"
    print(
        f"  anndata={'yes' if summary['anndata_available'] else 'no'}  "
        f"samtools={'yes' if summary['samtools_available'] else 'no'}  "
        f"mode={mode}  scope={scope}"
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
                             "budget. Default (unset) keeps strict mode. "
                             "All other classes stay strict regardless.")
    parser.add_argument(
        "--allow-subset",
        action="store_true",
        help="compare an explicitly curated, non-empty common subset (for the "
             "documented 1.x-to-2.0 validation). By default the five fixed "
             "pipeline_info outputs are required as proof of complete runs.",
    )
    args = parser.parse_args(argv)

    for d in (args.dir_a, args.dir_b):
        if not os.path.isdir(d):
            parser.error(f"not a directory: {d}")

    if args.envelope_max_flips is not None and args.envelope_max_flips < 0:
        parser.error("--envelope-max-flips must be >= 0")
    if args.max_diffs < 0:
        parser.error("--max-diffs must be >= 0")

    results, summary, failed = run(
        args.dir_a,
        args.dir_b,
        args.max_diffs,
        args.envelope_max_flips,
        require_complete=not args.allow_subset,
    )
    print_report(results, summary, args.max_diffs)

    if args.json_out:
        with open(args.json_out, "w", encoding="utf-8") as fh:
            json.dump({"summary": summary, "results": results}, fh, indent=2)
        print(f"JSON written to {args.json_out}")

    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
