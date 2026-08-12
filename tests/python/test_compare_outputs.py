"""Regression tests for the published-output equivalence comparator."""

import gzip
import importlib.util
import json
import re
import shutil
import subprocess
import sys
import tracemalloc
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parents[2]
COMPARATOR_PATH = REPO_ROOT / "tests" / "regression" / "compare_outputs.py"
SPEC = importlib.util.spec_from_file_location("compare_outputs", COMPARATOR_PATH)
compare_outputs = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(compare_outputs)


def _write(path, content):
    path.write_text(content, encoding="utf-8")
    return path


def _paired_paths(tmp_path, name):
    dir_a = tmp_path / "a"
    dir_b = tmp_path / "b"
    dir_a.mkdir(parents=True)
    dir_b.mkdir(parents=True)
    return dir_a / name, dir_b / name


def _write_gzip(path, content, *, mtime=0):
    with path.open("wb") as raw:
        with gzip.GzipFile(filename=path.name, fileobj=raw, mode="wb", mtime=mtime) as fh:
            fh.write(content.encode("utf-8"))
    return path


def _resolved_config(outdir, *, work_dir=None, genome="g"):
    config = {
        "input_csv": "/input/samples.csv",
        "outdir": outdir,
        "barcode_list_path": "/input/barcodes.csv",
        "barcode_correction_list_path": "/input/corrections.tsv",
        "barcode_pattern": "CCCCCCCCCCCCC",
        "minimum_count_threshold": 100,
        "sss_nmer": 8,
        "awsregion": None,
        "awsqueue": None,
        "max_memory": {"bytes": 274877906944},
        "max_cpus": 16,
        "max_time": {"durationInMillis": 864000000},
        "genome": genome,
    }
    if work_dir is not None:
        config["workDir"] = work_dir
    return config


def _write_required_pipeline_info(root):
    pipeline_info = root / "pipeline_info"
    pipeline_info.mkdir(parents=True, exist_ok=True)
    _write(
        pipeline_info / "resolved_configuration.txt",
        json.dumps(_resolved_config(str(root), work_dir=f"{root}/work")),
    )
    _write(
        pipeline_info / "execution_trace.txt",
        _trace_text(_trace_row(1, "GREET (sample-a)")),
    )
    _write(pipeline_info / "execution_report.html", _execution_report_html())
    _write(pipeline_info / "execution_timeline.html", _execution_timeline_html())
    _write(pipeline_info / "pipeline_dag.html", _pipeline_dag_html())


@pytest.mark.unit
@pytest.mark.parametrize(
    ("published_path", "expected_class"),
    [
        ("report/sample/sample.metrics.csv", "TEXT_EXACT"),
        ("report/multisample_out.csv", "TEXT_EXACT"),
        ("RSeQC/read_distribution/sample.raw_RSeQC.txt", "TEXT_EXACT"),
        ("deduplication/sample.dedup.log", "DEDUP_LOG"),
        ("pipeline_info/resolved_configuration.txt", "CONFIG"),
        ("pipeline_info/execution_trace.txt", "TRACE"),
        ("multiqc/sample_multiqc_data/multiqc.log", "MULTIQC_LOG"),
        ("multiqc/sample.multiqc.data.json", "MULTIQC_JSON"),
        ("multiqc/sample_multiqc_data/multiqc_data.json", "MULTIQC_JSON"),
        ("multiqc/sample_multiqc_data/multiqc_sources.txt", "MULTIQC_SOURCES"),
        ("count_matrix/raw/sample/barcodes.tsv.gz", "GZ_TEXT_EXACT"),
        ("count_matrix/filtered/sample/features.tsv.gz", "GZ_TEXT_EXACT"),
        ("count_matrix/raw/sample/matrix.mtx.gz", "MTX"),
        ("count_matrix/raw/sample/sample.raw_feature_bc_matrix.h5ad", "H5AD"),
        ("STAR/sample_Aligned.out.bam", "BAM"),
        ("featureCounts/sample.featureCounts.bam", "BAM"),
        ("deduplication/sample.dedup.bam", "BAM"),
        ("report/consolidated_report.html", "HTML"),
        ("report/multisample_qc_cascade.html", "HTML"),
        ("report/sample/sample.qc_cascade.html", "HTML"),
        ("plots/sample_pdf_with_cutoff.html", "HTML"),
        ("multiqc/multisample_multiqc.html", "HTML"),
        ("multiqc/single_sample_multiqc/sample/sample_multiqc.html", "HTML"),
        ("pipeline_info/execution_report.html", "HTML"),
        ("pipeline_info/execution_timeline.html", "HTML"),
        ("pipeline_info/pipeline_dag.html", "HTML"),
        ("multiqc/multisample_multiqc_data/multiqc_citations.txt", "BINARY_EXACT"),
        ("multiqc/multisample_multiqc_data/multiqc_rseqc.txt", "BINARY_EXACT"),
    ],
)
def test_classify_audits_the_configured_published_tree(published_path, expected_class):
    assert compare_outputs.classify(published_path) == expected_class


@pytest.mark.unit
def test_dispatch_is_defined_for_every_claimed_file_class():
    assert set(compare_outputs._COMPARATORS) == compare_outputs.EXACT_CLASSES


@pytest.mark.unit
def test_exact_text_and_binary_classes_fail_on_content_changes(tmp_path):
    text_root = tmp_path / "text"
    text_root.mkdir()
    text_a, text_b = _paired_paths(text_root, "sample.metrics.csv")
    _write(text_a, "metric,value\nreads,10\n")
    _write(text_b, "metric,value\nreads,11\n")
    verdict, detail = compare_outputs._compare_text_exact(text_a, text_b, 10)
    assert verdict == compare_outputs.DIFFER
    assert detail["diffs"][0]["line"] == 2

    text_a.write_bytes(b"\xff")
    text_b.write_bytes(b"\xff")
    verdict, detail = compare_outputs._compare_text_exact(text_a, text_b, 10)
    assert verdict == compare_outputs.DIFFER
    assert set(detail["validation_errors"]) == {"a", "b"}

    binary_root = tmp_path / "binary"
    binary_root.mkdir()
    binary_a, binary_b = _paired_paths(binary_root, "multiqc_citations.txt")
    binary_a.write_bytes(b"\x00stable-a")
    binary_b.write_bytes(b"\x00stable-b")
    verdict, detail = compare_outputs._compare_binary_exact(binary_a, binary_b, 10)
    assert verdict == compare_outputs.DIFFER
    assert detail["sha256_a"] != detail["sha256_b"]


def _dedup_log(input_reads=10, output_reads=8):
    return (
        f"INFO Reads: Input Reads: {input_reads}\n"
        f"INFO Number of reads out: {output_reads}\n"
    )


def _legacy_dedup_log(
    input_reads=10,
    output_reads=8,
    positions=8,
    date="2026-06-11",
    host="worker-a",
    pid=125,
    runtime=11,
    parsed_progress=None,
    written_progress=None,
    stable_overrides=None,
    extra_options=None,
):
    options = dict(compare_outputs._DEDUP_LEGACY_STABLE_OPTIONS)
    options.update({
        "stdin": "<_io.TextIOWrapper name='sample.bam' mode='r' encoding='UTF-8'>",
        "stdlog": (
            "<_io.TextIOWrapper name='sample.dedup.log' mode='a' "
            "encoding='UTF-8'>"
        ),
    })
    options.update(stable_overrides or {})
    options.update(extra_options or {})
    option_lines = [
        f"# {key:<40}: {value}" for key, value in sorted(options.items())
    ]
    if parsed_progress is None:
        parsed_progress = list(range(1_000_000, input_reads + 1, 1_000_000))
    if written_progress is None:
        written_progress = list(range(100_000, output_reads + 1, 100_000))
    info_lines = [
        f"{date} 15:34:20,558 INFO command: dedup --per-cell --in-sam "
        "-I sample.bam --log=sample.dedup.log",
        *(f"{date} 15:34:21,000 INFO Parsed {value} input reads"
          for value in parsed_progress),
        *(f"{date} 15:34:22,000 INFO Written out {value} reads"
          for value in written_progress),
        f"{date} 15:34:32,263 INFO Reads: Input Reads: {input_reads}",
        f"{date} 15:34:32,264 INFO Number of reads out: {output_reads}",
        f"{date} 15:34:32,264 INFO Total number of positions deduplicated: {positions}",
        f"{date} 15:34:32,264 INFO Mean number of unique UMIs per position: 1.00",
        f"{date} 15:34:32,264 INFO Max. number of unique UMIs per position: 1",
    ]
    return "\n".join([
        "# UMI-tools version: 1.1.2",
        "# output generated by dedup --per-cell --in-sam -I sample.bam --log=sample.dedup.log",
        f"# job started at Thu Jun 11 15:34:20 2026 on {host} -- uuid-a",
        f"# pid: {pid}, system: Linux test-kernel",
        *option_lines,
        *info_lines,
        f"# job finished in {runtime} seconds at Thu Jun 11 15:34:32 2026 -- 1.0 2.0 -- uuid-a",
        "",
    ])


@pytest.mark.unit
def test_dedup_log_compares_the_exact_current_summary(tmp_path):
    path_a, path_b = _paired_paths(tmp_path, "sample.dedup.log")
    _write(path_a, _dedup_log())
    _write(path_b, _dedup_log())
    assert compare_outputs._compare_dedup(path_a, path_b, 10)[0] == compare_outputs.EQUAL

    _write(path_b, _dedup_log(output_reads=7))
    assert compare_outputs._compare_dedup(path_a, path_b, 10)[0] == compare_outputs.DIFFER
    _write(path_b, "Input Reads: 10\n")
    verdict, detail = compare_outputs._compare_dedup(path_a, path_b, 10)
    assert verdict == compare_outputs.DIFFER
    assert "parse_error" in detail

    _write(path_b, "2026-08-11 pid=99\n" + _dedup_log())
    verdict, detail = compare_outputs._compare_dedup(path_a, path_b, 10)
    assert verdict == compare_outputs.DIFFER
    assert "validated completed" in detail["parse_error"]


@pytest.mark.unit
def test_dedup_log_validates_and_normalises_the_explicit_legacy_shape(tmp_path):
    path_a, path_b = _paired_paths(tmp_path, "sample.dedup.log")
    _write(path_a, _legacy_dedup_log())
    _write(
        path_b,
        _legacy_dedup_log(
            date="2026-06-12", host="worker-b", pid=999, runtime=27
        ),
    )
    verdict, detail = compare_outputs._compare_dedup(path_a, path_b, 10)
    assert verdict == compare_outputs.EQUAL
    assert detail == {"tuple": [10, 8]}

    # A validated legacy record and its current consolidated equivalent carry
    # the same stable read-count semantics for a 1.x -> 2.0 comparison.
    _write(path_b, _dedup_log())
    assert compare_outputs._compare_dedup(path_a, path_b, 10)[0] \
        == compare_outputs.EQUAL

    _write(path_b, _legacy_dedup_log(positions=7))
    verdict, detail = compare_outputs._compare_dedup(path_a, path_b, 10)
    assert verdict == compare_outputs.DIFFER
    assert "must equal reads out" in detail["parse_error"]


@pytest.mark.unit
def test_legacy_dedup_rejects_incoherent_progress_and_unknown_options(tmp_path):
    path_a, path_b = _paired_paths(tmp_path, "sample.dedup.log")
    _write(path_a, _legacy_dedup_log())

    _write(
        path_b,
        _legacy_dedup_log(parsed_progress=[999], written_progress=[777]),
    )
    verdict, detail = compare_outputs._compare_dedup(path_a, path_b, 10)
    assert verdict == compare_outputs.DIFFER
    assert "progress is inconsistent" in detail["parse_error"]

    _write(path_b, _legacy_dedup_log(extra_options={"invented": "value"}))
    verdict, detail = compare_outputs._compare_dedup(path_a, path_b, 10)
    assert verdict == compare_outputs.DIFFER
    assert "option set differs" in detail["parse_error"]

    _write(path_b, _legacy_dedup_log(stable_overrides={"threshold": "arbitrary"}))
    verdict, detail = compare_outputs._compare_dedup(path_a, path_b, 10)
    assert verdict == compare_outputs.DIFFER
    assert "threshold" in detail["parse_error"]


@pytest.mark.unit
def test_legacy_dedup_accepts_coherent_multi_record_progress(tmp_path):
    path_a, path_b = _paired_paths(tmp_path, "sample.dedup.log")
    content = _legacy_dedup_log(
        input_reads=2_000_001,
        output_reads=200_001,
        positions=200_001,
    )
    _write(path_a, content)
    _write(path_b, content)

    verdict, detail = compare_outputs._compare_dedup(path_a, path_b, 10)
    assert verdict == compare_outputs.EQUAL
    assert detail == {"tuple": [2_000_001, 200_001]}


@pytest.mark.unit
def test_dedup_log_accepts_current_normal_and_empty_producer_shapes(tmp_path):
    path_a, path_b = _paired_paths(tmp_path, "sample.dedup.log")
    normal = "INFO Reads: Input Reads: 10\nINFO Number of reads out: 8\n"
    _write(path_a, normal)
    _write(path_b, normal)
    verdict, detail = compare_outputs._compare_dedup(path_a, path_b, 10)
    assert verdict == compare_outputs.EQUAL
    assert detail == {"tuple": [10, 8]}

    empty = "INFO Reads: Input Reads: 0\nINFO Number of reads out: 0\n"
    _write(path_a, empty)
    _write(path_b, empty)
    verdict, detail = compare_outputs._compare_dedup(path_a, path_b, 10)
    assert verdict == compare_outputs.EQUAL
    assert detail == {"tuple": [0, 0]}

    # Pipeline 1.x's empty branch emitted one additional blank line. This exact
    # zero/zero legacy shape is accepted for upgrade validation, without
    # permitting generic trailing whitespace.
    legacy_empty = empty + "\n"
    _write(path_a, legacy_empty)
    _write(path_b, empty)
    verdict, detail = compare_outputs._compare_dedup(path_a, path_b, 10)
    assert verdict == compare_outputs.EQUAL
    assert detail == {"tuple": [0, 0]}

    _write(path_a, _dedup_log() + "\n")
    verdict, detail = compare_outputs._compare_dedup(path_a, path_b, 10)
    assert verdict == compare_outputs.DIFFER
    assert "validated completed" in detail["parse_error"]


@pytest.mark.unit
def test_dedup_log_rejects_impossible_or_duplicate_counts(tmp_path):
    path_a, path_b = _paired_paths(tmp_path, "sample.dedup.log")
    _write(path_a, _dedup_log())

    _write(path_b, _dedup_log(input_reads=5, output_reads=6))
    verdict, detail = compare_outputs._compare_dedup(path_a, path_b, 10)
    assert verdict == compare_outputs.DIFFER
    assert "cannot exceed" in detail["parse_error"]

    _write(path_b, _dedup_log() + "Input Reads: 10\n")
    verdict, detail = compare_outputs._compare_dedup(path_a, path_b, 10)
    assert verdict == compare_outputs.DIFFER
    assert "validated completed" in detail["parse_error"]


@pytest.mark.unit
def test_resolved_config_normalises_only_run_paths_and_fails_malformed_json(tmp_path):
    path_a, path_b = _paired_paths(tmp_path, "resolved_configuration.txt")
    _write(path_a, json.dumps(_resolved_config("s3://a", work_dir="/work/a")))
    # workDir is absent from genuine pinned Nextflow output unless configured.
    _write(path_b, json.dumps(_resolved_config("s3://b")))
    assert compare_outputs._compare_config(path_a, path_b, 10)[0] == compare_outputs.EQUAL

    _write(path_b, json.dumps(_resolved_config("s3://b", genome="other")))
    assert compare_outputs._compare_config(path_a, path_b, 10)[0] == compare_outputs.DIFFER
    _write(path_b, "not json")
    verdict, detail = compare_outputs._compare_config(path_a, path_b, 10)
    assert verdict == compare_outputs.DIFFER
    assert "parse_error" in detail

    for invalid in ("[]", '{"genome": NaN}'):
        _write(path_b, invalid)
        verdict, detail = compare_outputs._compare_config(path_a, path_b, 10)
        assert verdict == compare_outputs.DIFFER
        assert "parse_error" in detail


@pytest.mark.unit
@pytest.mark.parametrize(
    ("key", "value"),
    [
        ("outdir", 17),
        ("outdir", ""),
        ("workDir", None),
        ("workDir", ""),
        ("barcode_pattern", "not-a-pattern"),
        ("minimum_count_threshold", -1),
        ("max_cpus", 0),
        ("max_memory", {}),
        ("input_csv", 17),
    ],
)
def test_resolved_config_validates_volatile_fields_and_stable_payload(
    tmp_path, key, value
):
    valid, bad = _paired_paths(tmp_path, "resolved_configuration.txt")
    _write(valid, json.dumps(_resolved_config("s3://a")))
    invalid = _resolved_config("s3://b")
    invalid[key] = value
    _write(bad, json.dumps(invalid))

    verdict, detail = compare_outputs._compare_config(valid, bad, 10)

    assert verdict == compare_outputs.DIFFER
    assert "parse_error" in detail

    incomplete = _resolved_config("s3://b")
    incomplete.pop("barcode_pattern")
    _write(bad, json.dumps(incomplete))
    verdict, detail = compare_outputs._compare_config(valid, bad, 10)
    assert verdict == compare_outputs.DIFFER
    assert "stable pipeline parameters" in detail["parse_error"]


@pytest.mark.unit
def test_gz_text_compares_payload_not_gzip_headers_and_rejects_invalid_gzip(tmp_path):
    path_a, path_b = _paired_paths(tmp_path, "barcodes.tsv.gz")
    _write_gzip(path_a, "AAAC\nTTTG\n", mtime=1)
    _write_gzip(path_b, "AAAC\nTTTG\n", mtime=999)
    assert path_a.read_bytes() != path_b.read_bytes()
    assert compare_outputs._compare_gz_text_exact(path_a, path_b, 10)[0] == compare_outputs.EQUAL

    path_b.write_bytes(b"not gzip")
    verdict, detail = compare_outputs._compare_gz_text_exact(path_a, path_b, 10)
    assert verdict == compare_outputs.DIFFER
    assert "b" in detail["validation_errors"]


@pytest.mark.unit
def test_run_rejects_file_set_mismatch(tmp_path):
    dir_a = tmp_path / "a"
    dir_b = tmp_path / "b"
    dir_a.mkdir()
    dir_b.mkdir()
    _write_required_pipeline_info(dir_a)
    _write_required_pipeline_info(dir_b)
    _write(dir_a / "sample.metrics.csv", "metric,value\nreads,10\n")

    results, summary, failed = compare_outputs.run(dir_a, dir_b, 10)

    mismatches = [
        result for result in results
        if result["class"] == "FILE_SET_MISMATCH"
    ]
    assert mismatches == [{
        "path": "sample.metrics.csv",
        "class": "FILE_SET_MISMATCH",
        "verdict": compare_outputs.MISSING,
        "detail": {"present_in": "a", "required": False},
    }]
    assert summary["file_set_mismatch"] is True
    assert failed is True


@pytest.mark.unit
def test_cli_exit_codes_and_json_report(tmp_path):
    dir_a = tmp_path / "a"
    dir_b = tmp_path / "b"
    dir_a.mkdir()
    dir_b.mkdir()
    _write_required_pipeline_info(dir_a)
    _write_required_pipeline_info(dir_b)
    _write(dir_a / "sample.metrics.csv", "metric,value\nreads,10\n")
    _write(dir_b / "sample.metrics.csv", "metric,value\nreads,10\n")
    report = tmp_path / "comparison.json"

    equal = subprocess.run(
        [sys.executable, str(COMPARATOR_PATH), str(dir_a), str(dir_b), "--json", str(report)],
        check=False,
        capture_output=True,
        text=True,
    )
    assert equal.returncode == 0
    assert json.loads(report.read_text(encoding="utf-8"))["summary"]["failed"] is False

    _write(dir_b / "sample.metrics.csv", "metric,value\nreads,11\n")
    different = subprocess.run(
        [sys.executable, str(COMPARATOR_PATH), str(dir_a), str(dir_b)],
        check=False,
        capture_output=True,
        text=True,
    )
    assert different.returncode == 1
    assert "DIFFER" in different.stdout

    zero_detail = subprocess.run(
        [
            sys.executable, str(COMPARATOR_PATH), str(dir_a), str(dir_b),
            "--max-diffs", "0", "--json", str(report),
        ],
        check=False,
        capture_output=True,
        text=True,
    )
    assert zero_detail.returncode == 1
    metric_result = next(
        result for result in json.loads(report.read_text(encoding="utf-8"))["results"]
        if result["path"] == "sample.metrics.csv"
    )
    assert metric_result["detail"]["diffs"] == []

    negative_cap = subprocess.run(
        [
            sys.executable, str(COMPARATOR_PATH), str(dir_a), str(dir_b),
            "--max-diffs", "-1",
        ],
        check=False,
        capture_output=True,
        text=True,
    )
    assert negative_cap.returncode == 2
    assert "--max-diffs must be >= 0" in negative_cap.stderr

    # Curated cross-version comparisons are explicit. They do not silently
    # weaken the complete-tree default, and an empty subset still fails.
    subset_a = tmp_path / "subset-a"
    subset_b = tmp_path / "subset-b"
    subset_a.mkdir()
    subset_b.mkdir()
    _write(subset_a / "sample.metrics.csv", "metric,value\nreads,10\n")
    _write(subset_b / "sample.metrics.csv", "metric,value\nreads,10\n")
    subset = subprocess.run(
        [
            sys.executable, str(COMPARATOR_PATH), str(subset_a), str(subset_b),
            "--allow-subset",
        ],
        check=False,
        capture_output=True,
        text=True,
    )
    assert subset.returncode == 0
    assert "scope=explicit-subset" in subset.stdout

    empty_a = tmp_path / "empty-a"
    empty_b = tmp_path / "empty-b"
    empty_a.mkdir()
    empty_b.mkdir()
    empty = subprocess.run(
        [
            sys.executable, str(COMPARATOR_PATH), str(empty_a), str(empty_b),
            "--allow-subset",
        ],
        check=False,
        capture_output=True,
        text=True,
    )
    assert empty.returncode == 1
    assert "FILE_SET_MISMATCH" in empty.stdout


@pytest.mark.unit
def test_run_rejects_two_empty_or_identically_incomplete_output_trees(tmp_path):
    dir_a = tmp_path / "a"
    dir_b = tmp_path / "b"
    dir_a.mkdir()
    dir_b.mkdir()

    results, summary, failed = compare_outputs.run(dir_a, dir_b, 10)

    assert failed is True
    assert summary["file_set_mismatch"] is True
    assert {result["path"] for result in results} \
        == compare_outputs._REQUIRED_PUBLISHED_FILES
    assert all(result["detail"] == {
        "missing_from": ["a", "b"], "required": True,
    } for result in results)


def _trace_row(
    task_id,
    name,
    *,
    task_hash="ab/123456",
    native_id="1234",
    status="COMPLETED",
    exit_code=0,
    submit="2026-08-11 20:14:28.493",
    duration="121ms",
    realtime="4ms",
    cpu="62.7%",
    peak_rss="158.1 KB",
    peak_vmem="1 MB",
    rchar="2 KB",
    wchar="226 B",
):
    values = (
        task_id, task_hash, native_id, name, status, exit_code, submit,
        duration, realtime, cpu, peak_rss, peak_vmem, rchar, wchar,
    )
    return "\t".join(str(value) for value in values)


def _trace_text(*rows, header=None):
    columns = compare_outputs._TRACE_COLUMNS if header is None else header
    return "\t".join(columns) + "\n" + "\n".join(rows) + "\n"


@pytest.mark.unit
def test_trace_ignores_justified_volatility_and_row_order(tmp_path):
    path_a, path_b = _paired_paths(tmp_path, "execution_trace.txt")
    _write(
        path_a,
        _trace_text(
            _trace_row(1, "STAR (sample-a)"),
            _trace_row(2, "DEDUP (sample-a)", task_hash="cd/abcdef"),
        ),
    )
    _write(
        path_b,
        _trace_text(
            _trace_row(
                91,
                "DEDUP (sample-a)",
                task_hash="ef/654321",
                native_id="scheduler-b",
                submit="2026-08-12 01:00:00.000",
                duration="5s",
                realtime="4s",
                cpu="99%",
                peak_rss="2 GB",
                peak_vmem="3 GB",
                rchar="20 MB",
                wchar="10 MB",
            ),
            _trace_row(90, "STAR (sample-a)", task_hash="12/abcdef"),
        ),
    )

    verdict, detail = compare_outputs._compare_trace(path_a, path_b, 10)

    assert verdict == compare_outputs.EQUAL
    assert detail == {"task_count": 2}


@pytest.mark.unit
@pytest.mark.parametrize("difference", ["process_or_tag", "status", "exit", "multiplicity"])
def test_trace_rejects_stable_task_semantic_changes(tmp_path, difference):
    path_a, path_b = _paired_paths(tmp_path, "execution_trace.txt")
    _write(path_a, _trace_text(_trace_row(1, "STAR (sample-a)")))
    row_b = _trace_row(2, "STAR (sample-a)", task_hash="cd/abcdef")
    if difference == "process_or_tag":
        row_b = _trace_row(2, "STAR (sample-b)", task_hash="cd/abcdef")
    elif difference == "status":
        row_b = _trace_row(
            2, "STAR (sample-a)", task_hash="cd/abcdef", status="CACHED"
        )
    elif difference == "exit":
        row_b = _trace_row(
            2, "STAR (sample-a)", task_hash="cd/abcdef", exit_code=1
        )
    elif difference == "multiplicity":
        row_b += "\n" + _trace_row(3, "STAR (sample-a)", task_hash="ef/123456")
    _write(path_b, _trace_text(row_b))

    verdict, detail = compare_outputs._compare_trace(path_a, path_b, 10)

    assert verdict == compare_outputs.DIFFER
    assert "only_a" in detail or "execution_failures" in detail


@pytest.mark.unit
@pytest.mark.parametrize(
    "content",
    [
        "",
        "\t".join(compare_outputs._TRACE_COLUMNS) + "\n",
        _trace_text(_trace_row(1, "STAR (sample-a)"), header=("task_id", "name")),
        _trace_text("1\tab/123456\t123\tSTAR (sample-a)\tCOMPLETED\t0"),
        _trace_text(_trace_row("not-an-int", "STAR (sample-a)")),
        _trace_text(_trace_row(1, "STAR (sample-a)", task_hash="not-a-hash")),
        _trace_text(_trace_row(1, "STAR (sample-a)", status="RUNNING")),
        _trace_text(_trace_row(1, "STAR (sample-a)", duration="eventually")),
        _trace_text(
            _trace_row(1, "STAR (sample-a)"),
            _trace_row(1, "DEDUP (sample-a)", task_hash="cd/abcdef"),
        ),
    ],
)
def test_trace_fails_loud_on_malformed_or_incomplete_tsv(tmp_path, content):
    valid, invalid = _paired_paths(tmp_path, "execution_trace.txt")
    _write(valid, _trace_text(_trace_row(1, "STAR (sample-a)")))
    _write(invalid, content)

    verdict, detail = compare_outputs._compare_trace(valid, invalid, 10)

    assert verdict == compare_outputs.DIFFER
    assert "b" in detail["validation_errors"]


@pytest.mark.unit
def test_trace_rejects_unsuccessful_runs_even_when_both_match(tmp_path):
    path_a, path_b = _paired_paths(tmp_path, "execution_trace.txt")
    failed = _trace_text(
        _trace_row(1, "STAR (sample-a)", status="FAILED", exit_code=1)
    )
    _write(path_a, failed)
    _write(path_b, failed)

    verdict, detail = compare_outputs._compare_trace(path_a, path_b, 10)

    assert verdict == compare_outputs.DIFFER
    assert detail["execution_failures"]["a"][0]["exit"] == 1


def _multiqc_log(
    timestamp="2026-08-11 10:00:00,000",
    run_dir="/work/aa/111111",
    temp_dir="/tmp/tmp-a",
    found=1,
    include_update=False,
    level="INFO",
):
    def line(logger, record_level, message, millisecond="000"):
        current = timestamp[:-3] + millisecond
        return f"[{current}] {logger:<50} [{record_level:<7}]  {message}"

    records = [
        line("multiqc", "DEBUG", "This is MultiQC v1.14"),
        line(
            "multiqc",
            "DEBUG",
            "Command used: /usr/local/bin/multiqc . -f --title 'sample multiqc' "
            "--filename sample_multiqc.html -m unified_qc -m rseqc",
            "001",
        ),
    ]
    if include_update:
        records.append(
            line("multiqc", "WARNING", "MultiQC Version v1.35 now available!", "002")
        )
    records.extend([
        line("multiqc", "DEBUG", f"Working dir : {run_dir}", "003"),
        line(
            "multiqc",
            "DEBUG",
            f"Using temporary directory for creating report: {temp_dir}",
            "004",
        ),
        line("multiqc", "INFO", f"Search path : {run_dir}", "005"),
        line("multiqc.modules.rseqc", level, f"Found {found} reports", "006"),
        line("multiqc", "INFO", f"Report      : {run_dir}/sample_multiqc.html", "007"),
        line("multiqc", "INFO", f"Data        : {run_dir}/sample_multiqc_data", "008"),
        line(
            "multiqc",
            "DEBUG",
            f"Moving data file from '{temp_dir}/multiqc_data' "
            f"to '{run_dir}/sample_multiqc_data'",
            "009",
        ),
        line("multiqc", "INFO", "MultiQC complete", "010"),
    ])
    return "\n".join(records) + "\n"


@pytest.mark.unit
def test_multiqc_log_normalises_only_provenance_and_update_check(tmp_path):
    path_a, path_b = _paired_paths(tmp_path, "multiqc.log")
    _write(path_a, _multiqc_log())
    _write(
        path_b,
        _multiqc_log(
            timestamp="2026-08-12 11:22:33,000",
            run_dir="/work/bb/999999",
            temp_dir="/tmp/tmp-b",
            include_update=True,
        ),
    )

    verdict, detail = compare_outputs._compare_multiqc_log(path_a, path_b, 10)

    assert verdict == compare_outputs.EQUAL
    assert detail["record_count"] == 10


@pytest.mark.unit
def test_multiqc_log_rejects_stable_change_error_and_malformed_content(tmp_path):
    path_a, path_b = _paired_paths(tmp_path, "multiqc.log")
    _write(path_a, _multiqc_log())
    _write(path_b, _multiqc_log(found=2))
    assert compare_outputs._compare_multiqc_log(path_a, path_b, 10)[0] \
        == compare_outputs.DIFFER

    _write(path_b, _multiqc_log(level="ERROR"))
    verdict, detail = compare_outputs._compare_multiqc_log(path_a, path_b, 10)
    assert verdict == compare_outputs.DIFFER
    assert "logged ERROR" in detail["validation_errors"]["b"]

    _write(path_b, "arbitrary log text\n")
    verdict, detail = compare_outputs._compare_multiqc_log(path_a, path_b, 10)
    assert verdict == compare_outputs.DIFFER
    assert "malformed" in detail["validation_errors"]["b"]

    _write(path_b, _multiqc_log(run_dir=""))
    verdict, detail = compare_outputs._compare_multiqc_log(path_a, path_b, 10)
    assert verdict == compare_outputs.DIFFER
    assert "path must be non-empty" in detail["validation_errors"]["b"]


def _multiqc_json(run_dir="/work/aa/111111", creation="2026-08-11, 10:00 UTC"):
    return {
        "report_data_sources": {
            "RSeQC": {
                "read_distribution": {
                    "sample": f"{run_dir}/sample_rseqc_results.txt"
                }
            }
        },
        "report_general_stats_data": [{"sample": "sample", "reads": 10}],
        "report_general_stats_headers": {"reads": {"title": "Reads"}},
        "report_multiqc_command": (
            "/usr/local/bin/multiqc . -f --title 'sample multiqc' "
            "--filename sample_multiqc.html -m unified_qc -m rseqc"
        ),
        "report_plot_data": {"plot": {"datasets": [[1, 2, 3]]}},
        "report_saved_raw_data": {"rseqc": {"sample": {"reads": 10}}},
        "config_analysis_dir_abs": [run_dir],
        "config_analysis_dir": [run_dir],
        "config_creation_date": creation,
        "config_title": "sample multiqc",
        "config_version": "1.14",
    }


@pytest.mark.unit
def test_multiqc_json_normalises_paths_and_creation_time_only(tmp_path):
    path_a, path_b = _paired_paths(tmp_path, "multiqc_data.json")
    _write(path_a, json.dumps(_multiqc_json()))
    _write(
        path_b,
        json.dumps(_multiqc_json("/work/bb/999999", "2026-08-12, 11:22 UTC")),
    )

    assert compare_outputs._compare_multiqc_json(path_a, path_b, 10)[0] \
        == compare_outputs.EQUAL

    changed = _multiqc_json("/work/bb/999999", "2026-08-12, 11:22 UTC")
    changed["report_plot_data"]["plot"]["datasets"] = [[9, 9, 9]]
    _write(path_b, json.dumps(changed))
    verdict, detail = compare_outputs._compare_multiqc_json(path_a, path_b, 10)
    assert verdict == compare_outputs.DIFFER
    assert detail["differing_keys"] == ["report_plot_data"]


@pytest.mark.unit
@pytest.mark.parametrize(
    "invalid",
    [
        "not json",
        "{}",
        json.dumps({**_multiqc_json(), "config_version": "2.0"}),
        json.dumps({**_multiqc_json(), "config_version": "1.140"}),
        json.dumps({**_multiqc_json(), "report_plot_data": []}),
        json.dumps({**_multiqc_json(), "config_analysis_dir_abs": 17}),
        json.dumps({**_multiqc_json(), "config_analysis_dir": None}),
        json.dumps({**_multiqc_json(), "config_analysis_dir": [""]}),
        json.dumps({**_multiqc_json(), "report_data_sources": {
            "RSeQC": {"sample": ""}
        }}),
        json.dumps({**_multiqc_json(), "report_data_sources": {
            "RSeQC": {}
        }}),
    ],
)
def test_multiqc_json_fails_loud_on_invalid_contract(tmp_path, invalid):
    valid, bad = _paired_paths(tmp_path, "multiqc_data.json")
    _write(valid, json.dumps(_multiqc_json()))
    _write(bad, invalid)

    verdict, detail = compare_outputs._compare_multiqc_json(valid, bad, 10)

    assert verdict == compare_outputs.DIFFER
    assert "b" in detail["validation_errors"]


def _multiqc_sources(*rows):
    header = "\t".join(compare_outputs._MULTIQC_SOURCES_HEADER)
    return header + "\n" + "\n".join("\t".join(row) for row in rows) + "\n"


@pytest.mark.unit
def test_multiqc_sources_normalises_directories_and_row_order(tmp_path):
    path_a, path_b = _paired_paths(tmp_path, "multiqc_sources.txt")
    _write(
        path_a,
        _multiqc_sources(
            ("RSeQC", "read_distribution", "sample", "/work/a/sample_rseqc.txt"),
            ("fastp", "fastp", "sample", "/work/a/sample.fastp.json"),
        ),
    )
    _write(
        path_b,
        _multiqc_sources(
            ("fastp", "fastp", "sample", "/work/b/sample.fastp.json"),
            ("RSeQC", "read_distribution", "sample", "/work/b/sample_rseqc.txt"),
        ),
    )
    assert compare_outputs._compare_multiqc_sources(path_a, path_b, 10)[0] \
        == compare_outputs.EQUAL

    _write(
        path_b,
        _multiqc_sources(
            ("RSeQC", "read_distribution", "other", "/work/b/sample_rseqc.txt"),
            ("fastp", "fastp", "sample", "/work/b/sample.fastp.json"),
        ),
    )
    assert compare_outputs._compare_multiqc_sources(path_a, path_b, 10)[0] \
        == compare_outputs.DIFFER


@pytest.mark.unit
def test_multiqc_sources_rejects_malformed_tsv(tmp_path):
    valid, bad = _paired_paths(tmp_path, "multiqc_sources.txt")
    _write(
        valid,
        _multiqc_sources(
            ("RSeQC", "read_distribution", "sample", "/work/a/sample_rseqc.txt")
        ),
    )
    _write(bad, "Module\tSection\n")

    verdict, detail = compare_outputs._compare_multiqc_sources(valid, bad, 10)

    assert verdict == compare_outputs.DIFFER
    assert "b" in detail["validation_errors"]

    _write(bad, _multiqc_sources(("RSeQC", "read_distribution", "sample", "/")))
    verdict, detail = compare_outputs._compare_multiqc_sources(valid, bad, 10)
    assert verdict == compare_outputs.DIFFER
    assert "basename" in detail["validation_errors"]["b"]


def _plot_html(label="plot"):
    return (
        f"<div class='plotly-graph-div' id='{label}'></div>"
        f"<script>Plotly.newPlot('{label}', [], {{}})</script>"
    )


def _consolidated_html(run="run-a"):
    return (
        "<!doctype html><html><head>"
        f"<title>CS Genetics scRNA-seq report - {run}</title></head><body>"
        "<main id='main-content'><section id='overview'><h1>Overview</h1></section>"
        "<section id='cross-sample'><table class='cs-xsample-table'>"
        "<tr><th>sample</th></tr><tr><td>S1</td></tr></table></section></main>"
        "</body></html>"
    )


def _multiqc_html_config():
    """The semantic mqc_config shape emitted by pinned MultiQC 1.14."""
    return {
        "decimalPoint_format": None,
        "num_datasets_plot_limit": 50,
        "sample_names_rename": [],
        "show_hide_mode": [],
        "show_hide_patterns": [],
        "show_hide_regex": [],
        "thousandsSep_format": None,
    }


def _multiqc_html(config=None):
    config = _multiqc_html_config() if config is None else config
    body = config if isinstance(config, str) else json.dumps(config)
    return (
        "<!doctype html><html><head><title>Sample MultiQC Report</title></head>"
        "<body><script id='mqc_config' type='application/json'>"
        f"{body}</script><h1 id='page_title'>MultiQC</h1>"
        "<div id='mqc-module-section-rseqc'>RSeQC</div></body></html>"
    )


def _execution_report_html(
    task="GREET (sample-a)", status="COMPLETED", exit_code="0", summary_process=None
):
    process = task.partition(" (")[0]
    summary_process = process if summary_process is None else summary_process
    tag = task.rpartition(" (")[2][:-1] if task.endswith(")") and " (" in task else ""
    data = {
        "trace": [{
            "process": process,
            "name": task,
            "tag": tag,
            "status": status,
            "exit": exit_code,
            "task_id": "1",
            "workdir": "/volatile/work",
        }],
        "summary": [{"process": summary_process, "cpuUsage": {"mean": 10.0}}],
    }
    return (
        "<!doctype html><html><head><title>[run] Nextflow Workflow Report</title>"
        "</head><body><nav id='nf-report-navbar'>Report</nav>"
        "<table id='tasks_table'><tr><td>task</td></tr></table>"
        f"<script>window.data = {json.dumps(data)};</script></body></html>"
    )


def _execution_timeline_html(label="GREET (sample-a)", cached=False, offset=0):
    data = {
        "elapsed": "1.4s",
        "beginningMillis": 1000 + offset,
        "endingMillis": 2000 + offset,
        "processes": [{
            "label": label,
            "cached": cached,
            "index": 0,
            "times": [{"starting_time": 1100 + offset, "ending_time": 1900 + offset}],
        }],
    }
    return (
        "<!doctype html><html><body><h3>Processes execution timeline</h3>"
        "<span id='label_launch'>launch</span><span id='label_elapsed'>elapsed</span>"
        "<span id='label_legend'>legend</span><div id='timeline'>timeline</div>"
        f"<script>window.data = {json.dumps(data)}\n"
        ";\n\n"
        "  var prot = ((\"https:\" == document.location.protocol) "
        "? \"https://\" : \"http://\");\n"
        "  document.write(unescape(\"%3Clink href='\" + prot + "
        "\"fonts.googleapis.com/css?family=Lato' rel='stylesheet' "
        "type='text/css' %3E%3C/link%3E\"));\n"
        "</script></body></html>"
    )


def _pipeline_dag_html(process="GREET"):
    return (
        "<!doctype html><html><body><pre class='mermaid'>"
        "%%{init: {'theme': 'base'}}%%\nflowchart TB\n"
        f"v0([\"{process}\"])\nv0 --> v1\nv1[\"output\"]\n"
        "</pre><script type='module'>import mermaid from 'mermaid'; "
        "mermaid.initialize({startOnLoad: true});</script></body></html>"
    )


@pytest.mark.unit
def test_html_accepts_consolidated_report_parsed_elements(tmp_path):
    path_a, path_b = _paired_paths(tmp_path, "consolidated_report.html")
    _write(path_a, _consolidated_html("run-a"))
    _write(path_b, _consolidated_html("run-b"))

    verdict, _detail = compare_outputs._compare_html(path_a, path_b, 10)

    assert verdict == compare_outputs.PRESENT


@pytest.mark.unit
@pytest.mark.parametrize(
    ("filename", "left", "right"),
    [
        (
            "execution_report.html",
            _execution_report_html(),
            _execution_report_html(),
        ),
        (
            "execution_timeline.html",
            _execution_timeline_html(offset=0),
            _execution_timeline_html(offset=99_000),
        ),
        (
            "pipeline_dag.html",
            _pipeline_dag_html(),
            _pipeline_dag_html(),
        ),
    ],
)
def test_html_accepts_nextflow_26_04_1_structures_and_normalised_semantics(
    tmp_path, filename, left, right
):
    path_a, path_b = _paired_paths(tmp_path, filename)
    _write(path_a, left)
    _write(path_b, right)

    verdict, detail = compare_outputs._compare_html(path_a, path_b, 10)

    assert verdict == compare_outputs.EQUAL
    assert "Nextflow" in detail["reason"]


@pytest.mark.unit
@pytest.mark.parametrize(
    ("filename", "left", "right"),
    [
        (
            "execution_report.html",
            _execution_report_html("GREET (sample-a)"),
            _execution_report_html("OTHER (sample-a)"),
        ),
        (
            "execution_report.html",
            _execution_report_html("GREET (sample-a)"),
            _execution_report_html("GREET (sample-a)", summary_process="OTHER"),
        ),
        (
            "execution_timeline.html",
            _execution_timeline_html("GREET (sample-a)"),
            _execution_timeline_html("OTHER (sample-a)"),
        ),
        (
            "pipeline_dag.html",
            _pipeline_dag_html("GREET"),
            _pipeline_dag_html("OTHER"),
        ),
    ],
)
def test_html_rejects_nextflow_stable_semantic_changes(tmp_path, filename, left, right):
    path_a, path_b = _paired_paths(tmp_path, filename)
    _write(path_a, left)
    _write(path_b, right)

    verdict, detail = compare_outputs._compare_html(path_a, path_b, 10)

    assert verdict == compare_outputs.DIFFER
    assert detail["reason"] == "stable parsed Nextflow semantics differ"


@pytest.mark.unit
def test_execution_report_rejects_failed_task_data(tmp_path):
    path_a, path_b = _paired_paths(tmp_path, "execution_report.html")
    _write(path_a, _execution_report_html())
    _write(path_b, _execution_report_html(status="FAILED", exit_code="1"))

    verdict, detail = compare_outputs._compare_html(path_a, path_b, 10)

    assert verdict == compare_outputs.DIFFER
    assert "unsuccessful task" in detail["validation_errors"]["b"]


@pytest.mark.unit
@pytest.mark.parametrize(
    ("status", "exit_code"),
    [
        ("COMPLETED", 0.5),
        ("COMPLETED", -1),
        (["COMPLETED"], "0"),
    ],
)
def test_nextflow_report_rejects_malformed_stable_task_values(
    tmp_path, status, exit_code
):
    path_a, path_b = _paired_paths(tmp_path, "execution_report.html")
    _write(path_a, _execution_report_html())
    _write(path_b, _execution_report_html(status=status, exit_code=exit_code))

    verdict, detail = compare_outputs._compare_html(path_a, path_b, 10)

    assert verdict == compare_outputs.DIFFER
    assert "b" in detail["validation_errors"]


@pytest.mark.unit
@pytest.mark.parametrize(
    ("filename", "valid", "invalid"),
    [
        (
            "execution_report.html",
            _execution_report_html(),
            _execution_report_html().replace("id='tasks_table'", "id='not-tasks'"),
        ),
        (
            "execution_timeline.html",
            _execution_timeline_html(),
            _execution_timeline_html().replace(
                '"processes": [{', '"processes": [], "discarded": [{', 1
            ),
        ),
        (
            "pipeline_dag.html",
            _pipeline_dag_html(),
            (
                "<html><body><pre class='mermaid'>flowchart TB</pre>"
                "<script type='module'>import mermaid from 'mermaid';"
                "mermaid.initialize({});</script></body></html>"
            ),
        ),
    ],
)
def test_nextflow_html_fails_loud_on_missing_semantic_structure(
    tmp_path, filename, valid, invalid
):
    path_a, path_b = _paired_paths(tmp_path, filename)
    _write(path_a, valid)
    _write(path_b, invalid)

    verdict, detail = compare_outputs._compare_html(path_a, path_b, 10)

    assert verdict == compare_outputs.DIFFER
    assert "b" in detail["validation_errors"]


@pytest.mark.unit
@pytest.mark.parametrize(
    "filename",
    [
        "multisample_qc_cascade.html",
        "sample.qc_cascade.html",
        "sample_counts_pdf_with_threshold.html",
        "sample_barnyard_plot.html",
        "sample_pdf_with_cutoff.html",
        "sample_hsap_pdf_with_cutoff.html",
    ],
)
def test_html_audits_every_plot_filename_contract(tmp_path, filename):
    path_a, path_b = _paired_paths(tmp_path, filename)
    _write(path_a, _plot_html("a"))
    _write(path_b, _plot_html("b"))
    assert compare_outputs._compare_html(path_a, path_b, 10)[0] == compare_outputs.PRESENT


@pytest.mark.unit
def test_html_accepts_output_specific_plot_fragments(tmp_path):
    path_a, path_b = _paired_paths(tmp_path, "sample.qc_cascade.html")
    _write(path_a, _plot_html("volatile-a"))
    _write(path_b, _plot_html("volatile-b"))

    verdict, detail = compare_outputs._compare_html(path_a, path_b, 10)

    assert verdict == compare_outputs.PRESENT
    assert "output-specific" in detail["reason"]


@pytest.mark.unit
@pytest.mark.parametrize("filename", ["sample_multiqc.html", "multisample_multiqc.html"])
def test_html_accepts_multiqc_landmarks(tmp_path, filename):
    path_a, path_b = _paired_paths(tmp_path, filename)
    content = _multiqc_html()
    _write(path_a, content)
    _write(path_b, content)

    verdict, _detail = compare_outputs._compare_html(path_a, path_b, 10)

    assert verdict == compare_outputs.PRESENT


@pytest.mark.unit
def test_html_rejects_superficial_multiqc_marker_decoy(tmp_path):
    path_a, path_b = _paired_paths(tmp_path, "sample_multiqc.html")
    valid = _multiqc_html()
    _write(path_a, valid)
    _write(
        path_b,
        "<html><head><title>MultiQC Report</title></head>"
        "<body><div class='mqc_decoy'>arbitrary page</div></body></html>",
    )

    verdict, detail = compare_outputs._compare_html(path_a, path_b, 10)

    assert verdict == compare_outputs.DIFFER
    assert "configuration script" in detail["validation_errors"]["b"]


@pytest.mark.unit
@pytest.mark.parametrize(
    "invalid",
    [
        _multiqc_html("not json"),
        _multiqc_html({}),
        _multiqc_html({**_multiqc_html_config(), "num_datasets_plot_limit": None}),
        _multiqc_html({**_multiqc_html_config(), "show_hide_regex": "[]"}),
        _multiqc_html().replace("type='application/json'", "type='text/plain'"),
    ],
)
def test_multiqc_html_parses_and_validates_pinned_configuration(tmp_path, invalid):
    valid, bad = _paired_paths(tmp_path, "sample_multiqc.html")
    _write(valid, _multiqc_html())
    _write(bad, invalid)

    verdict, detail = compare_outputs._compare_html(valid, bad, 10)

    assert verdict == compare_outputs.DIFFER
    assert "MultiQC configuration" in detail["validation_errors"]["b"]


@pytest.mark.unit
@pytest.mark.parametrize(
    "content",
    [
        "",
        "   \n",
        "not an HTML document",
        "<div><span>broken</div>",
        "<div>",
        "<!doctype html><html><head></head><body></body></html>",
    ],
)
def test_html_rejects_empty_or_malformed_content(tmp_path, content):
    valid, invalid = _paired_paths(tmp_path, "sample.qc_cascade.html")
    _write(valid, _plot_html())
    _write(invalid, content)

    verdict, detail = compare_outputs._compare_html(valid, invalid, 10)

    assert verdict == compare_outputs.DIFFER
    assert "b" in detail["validation_errors"]


@pytest.mark.unit
def test_html_rejects_balanced_error_page_without_output_landmarks(tmp_path):
    path_a, path_b = _paired_paths(tmp_path, "sample.qc_cascade.html")
    _write(path_a, _plot_html())
    _write(
        path_b,
        "<!doctype html><html><body><h1>Report generation failed</h1></body></html>",
    )

    verdict, detail = compare_outputs._compare_html(path_a, path_b, 10)

    assert verdict == compare_outputs.DIFFER
    assert "error or failure page" in detail["validation_errors"]["b"]


@pytest.mark.unit
@pytest.mark.parametrize(
    "decoy",
    [
        "<!-- <div class='plotly-graph-div'></div>"
        "<script>Plotly.newPlot('x', [], {})</script> -->",
        "<p>&lt;div class='plotly-graph-div'&gt;&lt;/div&gt; "
        "Plotly.newPlot('x', [], {})</p>",
    ],
)
def test_html_ignores_comment_and_text_marker_decoys(tmp_path, decoy):
    path_a, path_b = _paired_paths(tmp_path, "sample.qc_cascade.html")
    _write(path_a, _plot_html())
    _write(path_b, f"<!doctype html><html><body>{decoy}</body></html>")

    verdict, detail = compare_outputs._compare_html(path_a, path_b, 10)

    assert verdict == compare_outputs.DIFFER
    assert "b" in detail["validation_errors"]


@pytest.mark.unit
@pytest.mark.parametrize(
    ("script_attrs", "script"),
    [
        ("", "/* Plotly.newPlot('plot', [], {}) */"),
        ("", "const marker = \"Plotly.newPlot('plot', [], {})\";"),
        ("", "Plotly.newPlot('different-element', [], {});"),
        ("", "if (true) /Plotly.newPlot('plot', [], {})/.test('x');"),
        ("", "if (true) {} /Plotly.newPlot('plot', [], {})/.test('x');"),
        ("", "while (true) { break\n/Plotly.newPlot('plot', [], {})/.test('x'); }"),
        (" type='application/json'", "Plotly.newPlot('plot', [], {});"),
        (" src='missing.js'", "Plotly.newPlot('plot', [], {});"),
        (" src", "Plotly.newPlot('plot', [], {});"),
    ],
)
def test_html_ignores_javascript_plotly_decoys_and_requires_target_link(
    tmp_path, script_attrs, script
):
    path_a, path_b = _paired_paths(tmp_path, "sample.qc_cascade.html")
    _write(path_a, _plot_html())
    _write(
        path_b,
        "<html><body><div class='plotly-graph-div' id='plot'></div>"
        f"<script{script_attrs}>{script}</script></body></html>",
    )

    verdict, detail = compare_outputs._compare_html(path_a, path_b, 10)

    assert verdict == compare_outputs.DIFFER
    assert "targeting a parsed Plotly graph" in detail["validation_errors"]["b"]


@pytest.mark.unit
@pytest.mark.parametrize("kind", ["report-comment", "report-string", "dag-comment", "dag-string"])
def test_nextflow_html_ignores_javascript_marker_decoys(tmp_path, kind):
    if kind.startswith("report"):
        filename = "execution_report.html"
        valid = _execution_report_html()
        data = json.dumps({
            "trace": [{
                "process": "GREET", "name": "GREET (sample-a)",
                "tag": "sample-a", "status": "COMPLETED", "exit": "0",
            }],
            "summary": [{"process": "GREET"}],
        })
        marker = f"window.data = {data};"
        script = f"/* {marker} */" if kind.endswith("comment") \
            else f"const decoy = {json.dumps(marker)};"
        invalid = re.sub(r"<script>.*?</script>", f"<script>{script}</script>", valid)
    else:
        filename = "pipeline_dag.html"
        valid = _pipeline_dag_html()
        marker = "import mermaid from 'mermaid'; mermaid.initialize({});"
        script = f"/* {marker} */" if kind.endswith("comment") \
            else f"const decoy = {json.dumps(marker)};"
        invalid = re.sub(
            r"<script type='module'>.*?</script>",
            f"<script type='module'>{script}</script>",
            valid,
        )
    path_a, path_b = _paired_paths(tmp_path, filename)
    _write(path_a, valid)
    _write(path_b, invalid)

    verdict, detail = compare_outputs._compare_html(path_a, path_b, 10)

    assert verdict == compare_outputs.DIFFER
    assert "b" in detail["validation_errors"]


@pytest.mark.unit
@pytest.mark.parametrize("script_attrs", [" type='application/json'", " src='missing.js'"])
def test_nextflow_window_data_requires_an_executable_inline_script(
    tmp_path, script_attrs
):
    path_a, path_b = _paired_paths(tmp_path, "execution_report.html")
    valid = _execution_report_html()
    invalid = valid.replace("<script>", f"<script{script_attrs}>")
    _write(path_a, valid)
    _write(path_b, invalid)

    verdict, detail = compare_outputs._compare_html(path_a, path_b, 10)

    assert verdict == compare_outputs.DIFFER
    assert "window.data" in detail["validation_errors"]["b"]


@pytest.mark.unit
def test_nextflow_window_data_rejects_regex_and_multiple_assignment_decoys(tmp_path):
    path_a, path_b = _paired_paths(tmp_path, "execution_report.html")
    valid = _execution_report_html()
    regex_decoy = re.sub(
        r"<script>.*?</script>",
        "<script>if (true) {} /window.data = "
        "{\"trace\":[],\"summary\":[]}/.test('x');</script>",
        valid,
    )
    _write(path_a, valid)
    _write(path_b, regex_decoy)
    verdict, detail = compare_outputs._compare_html(path_a, path_b, 10)
    assert verdict == compare_outputs.DIFFER
    assert "window.data" in detail["validation_errors"]["b"]

    assignment = re.search(r"<script>(.*?)</script>", valid).group(1)
    _write(
        path_b,
        valid.replace("</body>", f"<script>{assignment}</script></body>"),
    )
    verdict, detail = compare_outputs._compare_html(path_a, path_b, 10)
    assert verdict == compare_outputs.DIFFER
    assert "multiple executable" in detail["validation_errors"]["b"]

    _write(
        path_b,
        valid.replace("</body>", "<script>window.data = null;</script></body>"),
    )
    verdict, detail = compare_outputs._compare_html(path_a, path_b, 10)
    assert verdict == compare_outputs.DIFFER
    assert "must be an object" in detail["validation_errors"]["b"]


@pytest.mark.unit
@pytest.mark.parametrize(
    "replacement",
    [
        " && null;</script>",
        "; window.data.trace = [];</script>",
        " && null; window.data.trace = [];</script>",
    ],
)
def test_nextflow_window_data_rejects_expression_suffix_and_runtime_mutation(
    tmp_path, replacement
):
    path_a, path_b = _paired_paths(tmp_path, "execution_report.html")
    valid = _execution_report_html()
    invalid = valid.replace(";</script>", replacement)
    _write(path_a, valid)
    _write(path_b, invalid)

    verdict, detail = compare_outputs._compare_html(path_a, path_b, 10)

    assert verdict == compare_outputs.DIFFER
    assert "window.data" in detail["validation_errors"]["b"]


@pytest.mark.unit
@pytest.mark.parametrize(
    "executable_suffix",
    [
        '; window["data"] = null;',
        '; Object.defineProperty(window, "data", {value: null});',
        '; throw new Error("broken");',
    ],
)
def test_nextflow_window_data_rejects_every_unpinned_executable_trailer(
    tmp_path, executable_suffix
):
    path_a, path_b = _paired_paths(tmp_path, "execution_report.html")
    valid = _execution_report_html()
    invalid = valid.replace(";</script>", f"{executable_suffix}</script>")
    _write(path_a, valid)
    _write(path_b, invalid)

    verdict, detail = compare_outputs._compare_html(path_a, path_b, 10)

    assert verdict == compare_outputs.DIFFER
    assert "unexpected executable trailer" in detail["validation_errors"]["b"]


@pytest.mark.unit
def test_nextflow_window_data_trailers_are_filename_specific(tmp_path):
    report_a, report_b = _paired_paths(tmp_path / "report", "execution_report.html")
    valid_report = _execution_report_html()
    timeline_trailer = compare_outputs._NEXTFLOW_WINDOW_DATA_TRAILERS[
        "execution_timeline.html"
    ]
    _write(report_a, valid_report)
    _write(
        report_b,
        valid_report.replace(";</script>", f"{timeline_trailer}</script>"),
    )
    verdict, detail = compare_outputs._compare_html(report_a, report_b, 10)
    assert verdict == compare_outputs.DIFFER
    assert "unexpected executable trailer" in detail["validation_errors"]["b"]

    timeline_a, timeline_b = _paired_paths(
        tmp_path / "timeline", "execution_timeline.html"
    )
    valid_timeline = _execution_timeline_html()
    _write(timeline_a, valid_timeline)
    _write(timeline_b, valid_timeline.replace(timeline_trailer, ";"))
    verdict, detail = compare_outputs._compare_html(timeline_a, timeline_b, 10)
    assert verdict == compare_outputs.DIFFER
    assert "unexpected executable trailer" in detail["validation_errors"]["b"]


@pytest.mark.unit
def test_mermaid_initialiser_requires_an_executable_inline_module(tmp_path):
    path_a, path_b = _paired_paths(tmp_path, "pipeline_dag.html")
    valid = _pipeline_dag_html()
    invalid = valid.replace(
        "<script type='module'>",
        "<script type='module' src='missing.js'>",
    )
    _write(path_a, valid)
    _write(path_b, invalid)

    verdict, detail = compare_outputs._compare_html(path_a, path_b, 10)

    assert verdict == compare_outputs.DIFFER
    assert "Mermaid module" in detail["validation_errors"]["b"]


@pytest.mark.unit
def test_html_rejects_error_heading_even_with_valid_plot_elements(tmp_path):
    path_a, path_b = _paired_paths(tmp_path, "sample.qc_cascade.html")
    _write(path_a, _plot_html())
    _write(path_b, f"<html><body><h1>Internal Server Error</h1>{_plot_html()}</body></html>")

    verdict, detail = compare_outputs._compare_html(path_a, path_b, 10)

    assert verdict == compare_outputs.DIFFER
    assert "error or failure page" in detail["validation_errors"]["b"]


@pytest.mark.unit
def test_html_accepts_only_symmetric_named_empty_plot_sentinels(tmp_path):
    path_a, path_b = _paired_paths(tmp_path, "sample_pdf_with_cutoff.html")
    path_a.write_bytes(b"")
    path_b.write_bytes(b"")

    verdict, detail = compare_outputs._compare_html(path_a, path_b, 10)

    assert verdict == compare_outputs.PRESENT
    assert "sentinel" in detail["reason"]

    path_b.write_text(_plot_html(), encoding="utf-8")
    verdict, detail = compare_outputs._compare_html(path_a, path_b, 10)
    assert verdict == compare_outputs.DIFFER
    assert detail["validation_errors"]["a"] == "empty HTML"


@pytest.mark.unit
def test_invalid_html_drives_failed_run(tmp_path):
    dir_a = tmp_path / "a"
    dir_b = tmp_path / "b"
    dir_a.mkdir()
    dir_b.mkdir()
    _write_required_pipeline_info(dir_a)
    _write_required_pipeline_info(dir_b)
    _write(dir_a / "report.html", "")
    _write(dir_b / "report.html", "")

    results, summary, failed = compare_outputs.run(dir_a, dir_b, 10)

    report_result = next(result for result in results if result["path"] == "report.html")
    assert report_result["verdict"] == compare_outputs.DIFFER
    assert summary["failed"] is True
    assert failed is True


@pytest.mark.unit
def test_unrecognised_balanced_html_fails_closed(tmp_path):
    path_a, path_b = _paired_paths(tmp_path, "unexpected.html")
    content = "<!doctype html><html><body><main>looks fine</main></body></html>"
    _write(path_a, content)
    _write(path_b, content)

    verdict, detail = compare_outputs._compare_html(path_a, path_b, 10)

    assert verdict == compare_outputs.DIFFER
    assert set(detail["validation_errors"]) == {"a", "b"}


def _mtx_text(entries, *, shape=(2, 2), field="integer", banner=None, declared=None):
    if banner is None:
        banner = f"%%MatrixMarket matrix coordinate {field} general"
    if declared is None:
        declared = len(entries)
    lines = [banner, "% generated test matrix", f"{shape[0]} {shape[1]} {declared}"]
    lines.extend(entries)
    return "\n".join(lines) + "\n"


def _write_large_mtx(path, entry_count, *, prefix=(), field="integer"):
    """Write a sizeable valid archive without materialising all entries in RAM."""
    total = entry_count + len(prefix)
    with path.open("wb") as raw:
        with gzip.GzipFile(
            filename=path.name, fileobj=raw, mode="wb", mtime=0
        ) as compressed:
            compressed.write(
                (
                    f"%%MatrixMarket matrix coordinate {field} general\n"
                    f"4096 64 {total}\n"
                    + "".join(f"{entry}\n" for entry in prefix)
                ).encode("ascii")
            )
            for start in range(0, entry_count, 10_000):
                stop = min(entry_count, start + 10_000)
                block = "".join(
                    f"{10 + (index % 4087)} "
                    f"{1 + ((index // 4087) % 64)} "
                    f"{1 + (index % 7)}\n"
                    for index in range(start, stop)
                )
                compressed.write(block.encode("ascii"))
    return path


@pytest.mark.unit
@pytest.mark.parametrize(
    "archive_kind", ["barcodes", "features", "matrix-header", "matrix-order"]
)
def test_complete_run_gates_tripartite_archive_bytes_after_validation(
    tmp_path, archive_kind
):
    dir_a = tmp_path / "a"
    dir_b = tmp_path / "b"
    dir_a.mkdir()
    dir_b.mkdir()
    _write_required_pipeline_info(dir_a)
    _write_required_pipeline_info(dir_b)
    relative_dir = Path("count_matrix/raw_feature_bc_matrix/SAMPLE")
    (dir_a / relative_dir).mkdir(parents=True)
    (dir_b / relative_dir).mkdir(parents=True)

    if archive_kind in {"barcodes", "features"}:
        relative = relative_dir / f"{archive_kind}.tsv.gz"
        payload = "AAAC\nTTTG\n" if archive_kind == "barcodes" \
            else "gene-1\tGene 1\n"
        _write_gzip(dir_a / relative, payload, mtime=1)
        _write_gzip(dir_b / relative, payload, mtime=999)
    else:
        relative = relative_dir / "matrix.mtx.gz"
        entries_a = ["1 1 2", "2 1 3", "2 2 4"]
        entries_b = entries_a if archive_kind == "matrix-header" \
            else list(reversed(entries_a))
        _write_gzip(
            dir_a / relative,
            _mtx_text(entries_a),
            mtime=1,
        )
        _write_gzip(
            dir_b / relative,
            _mtx_text(entries_b),
            mtime=999,
        )

    results, _summary, failed = compare_outputs.run(dir_a, dir_b, 10)
    archive_result = next(result for result in results if result["path"] == str(relative))
    assert failed is True
    assert archive_result["verdict"] == compare_outputs.DIFFER
    assert "gzip archive bytes differ" in archive_result["detail"]["reason"]

    # Curated 1.x-to-2.0 subset diagnostics validate the legacy decompressed
    # payload contract, but do not claim deterministic legacy gzip headers.
    subset_results, _summary, subset_failed = compare_outputs.run(
        dir_a, dir_b, 10, require_complete=False
    )
    subset_archive = next(
        result for result in subset_results if result["path"] == str(relative)
    )
    assert subset_failed is False
    assert subset_archive["verdict"] == compare_outputs.EQUAL


@pytest.mark.unit
def test_mtx_is_order_insensitive_after_strict_validation(tmp_path):
    path_a, path_b = _paired_paths(tmp_path, "matrix.mtx.gz")
    _write_gzip(path_a, _mtx_text(["1 1 2", "2 1 3", "2 2 4"]), mtime=1)
    _write_gzip(path_b, _mtx_text(["2 2 4", "1 1 2", "2 1 3"]), mtime=999)

    verdict, detail = compare_outputs._compare_mtx(path_a, path_b, 10)

    assert verdict == compare_outputs.EQUAL
    assert detail == {}


@pytest.mark.unit
def test_mtx_real_values_use_exact_decimal_numeric_semantics(tmp_path):
    path_a, path_b = _paired_paths(tmp_path, "matrix.mtx.gz")
    _write_gzip(
        path_a,
        _mtx_text(["1 1 1", "2 1 +01.000e+0"], field="real"),
    )
    _write_gzip(
        path_b,
        _mtx_text(["2 1 1e0", "1 1 1.0"], field="real"),
    )
    assert compare_outputs._compare_mtx(path_a, path_b, 10)[0] \
        == compare_outputs.EQUAL

    _write_gzip(
        path_a,
        _mtx_text(["1 1 9007199254740992"], field="real", declared=1),
    )
    _write_gzip(
        path_b,
        _mtx_text(["1 1 9007199254740993"], field="real", declared=1),
    )
    assert compare_outputs._compare_mtx(path_a, path_b, 10)[0] \
        == compare_outputs.DIFFER


@pytest.mark.unit
def test_mtx_real_rejects_underflowed_negative_and_python_number_syntax(tmp_path):
    valid, bad = _paired_paths(tmp_path, "matrix.mtx.gz")
    _write_gzip(valid, _mtx_text(["1 1 0"], field="real", declared=1))
    for token, message in (("-1e-4000", "non-negative"), ("1_0", "invalid real")):
        _write_gzip(bad, _mtx_text([f"1 1 {token}"], field="real", declared=1))
        verdict, detail = compare_outputs._compare_mtx(valid, bad, 10)
        assert verdict == compare_outputs.DIFFER
        assert message in detail["validation_errors"]["b"]


@pytest.mark.unit
def test_mtx_envelope_accepts_only_net_preserving_support_changes(tmp_path):
    path_a, path_b = _paired_paths(tmp_path, "matrix.mtx.gz")
    _write_gzip(path_a, _mtx_text(["1 1 1", "2 1 1"]))
    _write_gzip(path_b, _mtx_text(["2 1 2"], declared=1))

    verdict, detail = compare_outputs._compare_mtx(
        path_a, path_b, 10, envelope_max_flips=2
    )

    # A MatrixMarket archive's declared nnz is structural, but support changes
    # are still valid envelope candidates just like H5AD support changes.
    assert verdict == compare_outputs.ENVELOPE_OK
    assert detail["col_sums_preserved"] is True
    assert detail["n_diff_entries"] == 2

    _write_gzip(path_b, _mtx_text(["2 1 3"], declared=1))
    verdict, detail = compare_outputs._compare_mtx(
        path_a, path_b, 10, envelope_max_flips=99
    )
    assert verdict == compare_outputs.DIFFER
    assert detail["col_sums_preserved"] is False

    _write_gzip(
        path_a,
        _mtx_text(["1 1 0.1", "2 1 0.2"], field="real"),
    )
    _write_gzip(
        path_b,
        _mtx_text(["2 1 0.3"], field="real", declared=1),
    )
    verdict, detail = compare_outputs._compare_mtx(
        path_a, path_b, 10, envelope_max_flips=2
    )
    assert verdict == compare_outputs.ENVELOPE_OK
    assert detail["col_sums_preserved"] is True


@pytest.mark.unit
@pytest.mark.parametrize(
    ("invalid_text", "expected"),
    [
        ("", "empty MatrixMarket"),
        ("%MatrixMarket matrix coordinate integer general\n1 1 0\n", "banner"),
        (_mtx_text([], banner="%%MatrixMarket vector coordinate integer general"), "object"),
        (_mtx_text([], banner="%%MatrixMarket matrix array integer general"), "format"),
        (_mtx_text([], banner="%%MatrixMarket matrix coordinate complex general"), "field"),
        (_mtx_text([], banner="%%MatrixMarket matrix coordinate integer symmetric"), "symmetry"),
        ("%%MatrixMarket matrix coordinate integer general\n2 2\n", "dimensions"),
        (_mtx_text(["1 1"], declared=1), "coordinate entry"),
        (_mtx_text(["0 1 1"], declared=1), "1-based bounds"),
        (_mtx_text(["3 1 1"], declared=1), "1-based bounds"),
        (_mtx_text(["1 1 1.5"], declared=1), "integer value"),
        (_mtx_text(["1 1 -1"], declared=1), "non-negative"),
        (_mtx_text(["1 1 1"], declared=2), "declared nnz"),
        (_mtx_text(["1 1 nan"], field="real", declared=1), "non-finite"),
        (_mtx_text(["1 1 inf"], field="real", declared=1), "non-finite"),
        (_mtx_text(["1 1 1_0"], field="real", declared=1), "invalid real"),
        (_mtx_text(["1 1 1e10001"], field="real", declared=1), "safety limit"),
        (
            _mtx_text(["1 1 1e999999999", "2 1 1"], field="real"),
            "safety limit",
        ),
    ],
)
def test_mtx_fails_loud_on_invalid_matrixmarket_contract(tmp_path, invalid_text, expected):
    valid, invalid = _paired_paths(tmp_path, "matrix.mtx.gz")
    _write_gzip(valid, _mtx_text(["1 1 1"], declared=1))
    _write_gzip(invalid, invalid_text)

    verdict, detail = compare_outputs._compare_mtx(valid, invalid, 10)

    assert verdict == compare_outputs.DIFFER
    assert expected in detail["validation_errors"]["b"]


@pytest.mark.unit
@pytest.mark.parametrize("location", ["dimension", "row", "value"])
def test_mtx_huge_integer_tokens_fail_as_validation_differences(tmp_path, location):
    valid, bad = _paired_paths(tmp_path, "matrix.mtx.gz")
    huge = "9" * (compare_outputs._MTX_MAX_INTEGER_DIGITS + 1)
    _write_gzip(valid, _mtx_text(["1 1 1"], declared=1))
    if location == "dimension":
        invalid = (
            "%%MatrixMarket matrix coordinate integer general\n"
            f"{huge} 1 0\n"
        )
    elif location == "row":
        invalid = _mtx_text([f"{huge} 1 1"], declared=1)
    else:
        invalid = _mtx_text([f"1 1 {huge}"], declared=1)
    _write_gzip(bad, invalid)

    verdict, detail = compare_outputs._compare_mtx(valid, bad, 10)

    assert verdict == compare_outputs.DIFFER
    assert "safety limit" in detail["validation_errors"]["b"]


@pytest.mark.unit
@pytest.mark.parametrize("invalid_kind", ["arbitrary", "truncated"])
def test_mtx_rejects_corrupt_or_truncated_gzip(tmp_path, invalid_kind):
    valid, invalid = _paired_paths(tmp_path, "matrix.mtx.gz")
    _write_gzip(valid, _mtx_text(["1 1 1"], declared=1))
    if invalid_kind == "arbitrary":
        invalid.write_bytes(b"not gzip")
    else:
        invalid.write_bytes(valid.read_bytes()[:-5])

    verdict, detail = compare_outputs._compare_mtx(valid, invalid, 10)

    assert verdict == compare_outputs.DIFFER
    assert "b" in detail["validation_errors"]


@pytest.mark.unit
def test_mtx_envelope_sums_large_finite_reals_exactly_without_float_overflow(tmp_path):
    path_a, path_b = _paired_paths(tmp_path, "matrix.mtx.gz")
    _write_gzip(
        path_a,
        _mtx_text(["1 1 1e308", "2 1 1e308"], shape=(2, 1), field="real"),
    )
    _write_gzip(
        path_b,
        _mtx_text(["1 1 1.1e308", "2 1 9e307"], shape=(2, 1), field="real"),
    )

    verdict, detail = compare_outputs._compare_mtx(
        path_a, path_b, 10, envelope_max_flips=10
    )

    assert verdict == compare_outputs.ENVELOPE_OK
    assert detail["col_sums_preserved"] is True


@pytest.mark.unit
def test_mtx_forced_external_sort_preserves_multiplicity_envelope_and_fan_in(
    tmp_path, monkeypatch
):
    path_a, path_b = _paired_paths(tmp_path, "matrix.mtx.gz")
    monkeypatch.setattr(compare_outputs, "_MTX_SORT_CHUNK_BYTES", 16 * 1024)
    monkeypatch.setattr(compare_outputs, "_MTX_MERGE_FAN_IN", 3)
    merge_group_sizes = []
    final_run_counts = []
    original_merge = compare_outputs._merge_mtx_run_group
    original_iter = compare_outputs._iter_mtx_records

    def guarded_merge(paths, output):
        merge_group_sizes.append(len(paths))
        assert len(paths) <= 3
        return original_merge(paths, output)

    def guarded_iter(matrix):
        final_run_counts.append(len(matrix["runs"]))
        assert len(matrix["runs"]) <= 3
        yield from original_iter(matrix)

    monkeypatch.setattr(compare_outputs, "_merge_mtx_run_group", guarded_merge)
    monkeypatch.setattr(compare_outputs, "_iter_mtx_records", guarded_iter)

    _write_large_mtx(path_a, 30_000, prefix=("1 1 1", "1 1 1"))
    _write_large_mtx(path_b, 30_000, prefix=("1 1 1",))
    verdict, detail = compare_outputs._compare_mtx(path_a, path_b, 10)
    assert verdict == compare_outputs.DIFFER
    assert detail["n_only_a"] == 1

    _write_large_mtx(path_a, 30_000, prefix=("1 1 1", "2 1 1"))
    _write_large_mtx(path_b, 30_000, prefix=("2 1 2",))
    verdict, detail = compare_outputs._compare_mtx(
        path_a, path_b, 10, envelope_max_flips=2
    )
    assert verdict == compare_outputs.ENVELOPE_OK
    assert detail["n_diff_entries"] == 2
    assert detail["col_sums_preserved"] is True
    assert merge_group_sizes and final_run_counts


@pytest.mark.unit
def test_mtx_paired_large_comparison_memory_is_scale_bounded(tmp_path, monkeypatch):
    monkeypatch.setattr(compare_outputs, "_MTX_SORT_CHUNK_BYTES", 64 * 1024)
    monkeypatch.setattr(compare_outputs, "_MTX_MERGE_FAN_IN", 4)
    small_a, small_b = _paired_paths(tmp_path / "small", "matrix.mtx.gz")
    large_a, large_b = _paired_paths(tmp_path / "large", "matrix.mtx.gz")
    _write_large_mtx(small_a, 40_000)
    _write_large_mtx(small_b, 40_000)
    _write_large_mtx(large_a, 200_000)
    _write_large_mtx(large_b, 200_000)

    def peak_for(left, right):
        tracemalloc.start()
        try:
            verdict, _detail = compare_outputs._compare_mtx(left, right, 1)
            _current, peak = tracemalloc.get_traced_memory()
        finally:
            tracemalloc.stop()
        assert verdict == compare_outputs.EQUAL
        return peak

    small_peak = peak_for(small_a, small_b)
    large_peak = peak_for(large_a, large_b)

    # Five times as many entries must not cause entry-proportional Python RAM.
    assert large_peak < small_peak * 2
    assert large_peak < 16 * 1024 * 1024


@pytest.mark.unit
def test_mtx_external_sort_scratch_is_removed_after_parse_error(tmp_path, monkeypatch):
    valid, bad = _paired_paths(tmp_path, "matrix.mtx.gz")
    _write_large_mtx(valid, 30_000)
    _write_large_mtx(bad, 30_000)
    # Leave a syntactically valid gzip stream whose declared count is wrong.
    with gzip.open(bad, "rt", encoding="ascii") as fh:
        invalid_text = fh.read().replace("4096 64 30000", "4096 64 30001", 1)
    _write_gzip(bad, invalid_text)
    monkeypatch.setattr(compare_outputs, "_MTX_SORT_CHUNK_BYTES", 16 * 1024)

    scratch_parent = tmp_path / "comparator-scratch"
    scratch_parent.mkdir()
    real_temporary_directory = compare_outputs.tempfile.TemporaryDirectory

    class TrackedTemporaryDirectory(real_temporary_directory):
        def __init__(self, *args, **kwargs):
            kwargs["dir"] = scratch_parent
            super().__init__(*args, **kwargs)

    monkeypatch.setattr(
        compare_outputs.tempfile,
        "TemporaryDirectory",
        TrackedTemporaryDirectory,
    )
    verdict, detail = compare_outputs._compare_mtx(valid, bad, 10)

    assert verdict == compare_outputs.DIFFER
    assert "declared nnz" in detail["validation_errors"]["b"]
    assert list(scratch_parent.iterdir()) == []


def _make_adata(x, *, called=(True, False)):
    anndata = pytest.importorskip("anndata")
    np = pytest.importorskip("numpy")
    pd = pytest.importorskip("pandas")
    sp = pytest.importorskip("scipy.sparse")

    obs = pd.DataFrame(
        {"is_single_cell": list(called), "total_counts": np.asarray(x).sum(axis=1)},
        index=["sample_AAAC", "sample_TTTG"],
    )
    var = pd.DataFrame(
        {"gene_id": ["gene-1", "gene:2"], "is_mito": [False, True]},
        index=["Gene One", "Gene Two"],
    )
    return anndata.AnnData(sp.csr_matrix(np.asarray(x, dtype=np.float32)), obs, var)


def _make_x_only_adata(x, *, dtype):
    anndata = pytest.importorskip("anndata")
    np = pytest.importorskip("numpy")
    pd = pytest.importorskip("pandas")
    sp = pytest.importorskip("scipy.sparse")

    values = np.asarray(x, dtype=dtype)
    obs = pd.DataFrame(
        {"total_counts": np.zeros(values.shape[0], dtype=np.float64)},
        index=[f"barcode-{index}" for index in range(values.shape[0])],
    )
    var = pd.DataFrame(index=[f"gene-{index}" for index in range(values.shape[1])])
    return anndata.AnnData(sp.csr_matrix(values), obs=obs, var=var)


@pytest.mark.integration
def test_h5ad_compares_loaded_anndata_semantics(tmp_path):
    path_a = tmp_path / "a.h5ad"
    path_b = tmp_path / "b.h5ad"
    _make_adata([[1, 0], [2, 3]]).write_h5ad(path_a, compression="gzip")
    _make_adata([[1, 0], [2, 3]]).write_h5ad(path_b)

    verdict, detail = compare_outputs._compare_h5ad(path_a, path_b, 10)

    assert verdict == compare_outputs.EQUAL
    assert detail == {"shape": [2, 2], "nnz": 3}


@pytest.mark.integration
def test_h5ad_matrix_comparison_is_sparse_dense_storage_neutral(tmp_path):
    path_a = tmp_path / "a.h5ad"
    path_b = tmp_path / "b.h5ad"
    adata_a = _make_adata([[1, 0], [2, 3]])
    adata_b = _make_adata([[1, 0], [2, 3]])
    adata_b.X = adata_b.X.toarray()
    adata_a.write_h5ad(path_a)
    adata_b.write_h5ad(path_b)

    verdict, detail = compare_outputs._compare_h5ad(path_a, path_b, 10)

    assert verdict == compare_outputs.EQUAL
    assert detail == {"shape": [2, 2], "nnz": 3}


@pytest.mark.integration
@pytest.mark.parametrize("difference", ["matrix", "annotation", "layer"])
def test_h5ad_rejects_semantic_differences(tmp_path, difference):
    path_a = tmp_path / "a.h5ad"
    path_b = tmp_path / "b.h5ad"
    adata_a = _make_adata([[1, 0], [2, 3]])
    adata_b = _make_adata([[1, 0], [2, 3]])
    if difference == "matrix":
        adata_b.X[0, 0] = 9
    elif difference == "annotation":
        adata_b.obs["is_single_cell"] = [False, False]
    else:
        adata_b.layers["counts"] = adata_b.X.copy()
    adata_a.write_h5ad(path_a)
    adata_b.write_h5ad(path_b)

    verdict, detail = compare_outputs._compare_h5ad(path_a, path_b, 10)

    assert verdict == compare_outputs.DIFFER
    assert "mismatches" in detail


@pytest.mark.integration
def test_h5ad_envelope_still_allows_net_preserving_gene_flip(tmp_path):
    path_a = tmp_path / "a.h5ad"
    path_b = tmp_path / "b.h5ad"
    _make_adata([[1, 0], [2, 3]]).write_h5ad(path_a)
    _make_adata([[0, 1], [2, 3]]).write_h5ad(path_b)

    verdict, detail = compare_outputs._compare_h5ad(
        path_a, path_b, 10, envelope_max_flips=1
    )

    assert verdict == compare_outputs.ENVELOPE_OK
    assert detail["col_sums_preserved"] is True
    assert detail["n_diff_entries"] == 1


@pytest.mark.integration
def test_h5ad_envelope_allows_net_preserving_support_change(tmp_path):
    path_a = tmp_path / "a.h5ad"
    path_b = tmp_path / "b.h5ad"
    _make_adata([[1, 1], [2, 3]]).write_h5ad(path_a)
    _make_adata([[0, 2], [2, 3]]).write_h5ad(path_b)

    verdict, detail = compare_outputs._compare_h5ad(
        path_a, path_b, 10, envelope_max_flips=2
    )

    assert verdict == compare_outputs.ENVELOPE_OK
    assert detail["col_sums_preserved"] is True
    assert detail["n_diff_entries"] == 2


@pytest.mark.integration
@pytest.mark.parametrize("nonfinite", [float("nan"), float("inf"), float("-inf")])
def test_h5ad_rejects_nonfinite_x_even_when_both_files_match(tmp_path, nonfinite):
    np = pytest.importorskip("numpy")
    path_a = tmp_path / "a.h5ad"
    path_b = tmp_path / "b.h5ad"
    _make_x_only_adata([[nonfinite, 1.0, 0.0]], dtype=np.float64).write_h5ad(path_a)
    _make_x_only_adata([[nonfinite, 1.0, 0.0]], dtype=np.float64).write_h5ad(path_b)

    verdict, detail = compare_outputs._compare_h5ad(path_a, path_b, 10)

    assert verdict == compare_outputs.DIFFER
    assert set(detail["validation_errors"]) == {"a", "b"}
    assert all("non-finite" in error for error in detail["validation_errors"].values())


@pytest.mark.integration
def test_h5ad_nonfinite_cannot_spoof_net_preserving_envelope(tmp_path):
    np = pytest.importorskip("numpy")
    path_a = tmp_path / "a.h5ad"
    path_b = tmp_path / "b.h5ad"
    _make_x_only_adata([[np.inf, 1.0, 0.0]], dtype=np.float64).write_h5ad(path_a)
    _make_x_only_adata([[np.inf, 0.0, 2.0]], dtype=np.float64).write_h5ad(path_b)

    verdict, detail = compare_outputs._compare_h5ad(
        path_a, path_b, 10, envelope_max_flips=99
    )

    assert verdict == compare_outputs.DIFFER
    assert set(detail["validation_errors"]) == {"a", "b"}


@pytest.mark.integration
@pytest.mark.parametrize("location", ["sparse-X", "dense-X", "raw-X"])
def test_h5ad_rejects_negative_count_values_in_every_x_path(tmp_path, location):
    np = pytest.importorskip("numpy")
    sp = pytest.importorskip("scipy.sparse")
    path_a = tmp_path / "a.h5ad"
    path_b = tmp_path / "b.h5ad"

    for path in (path_a, path_b):
        adata = _make_x_only_adata([[-1.0, 2.0]], dtype=np.float64)
        if location == "dense-X":
            adata.X = adata.X.toarray()
        elif location == "raw-X":
            adata.raw = adata.copy()
            adata.X = sp.csr_matrix(np.asarray([[1.0, 2.0]]))
        adata.write_h5ad(path)

    verdict, detail = compare_outputs._compare_h5ad(
        path_a, path_b, 10, envelope_max_flips=99
    )

    assert verdict == compare_outputs.DIFFER
    assert set(detail["validation_errors"]) == {"a", "b"}
    assert all("negative count" in error
               for error in detail["validation_errors"].values())


@pytest.mark.integration
def test_h5ad_envelope_rejects_float_summation_overflow(tmp_path):
    np = pytest.importorskip("numpy")
    path_a = tmp_path / "a.h5ad"
    path_b = tmp_path / "b.h5ad"
    _make_x_only_adata([[1e308, 1e308]], dtype=np.float64).write_h5ad(path_a)
    _make_x_only_adata([[1.1e308, 9e307]], dtype=np.float64).write_h5ad(path_b)

    verdict, detail = compare_outputs._compare_h5ad(
        path_a, path_b, 10, envelope_max_flips=99
    )

    assert verdict == compare_outputs.DIFFER
    assert "overflow" in detail["validation_errors"]["envelope"]


@pytest.mark.integration
def test_sparse_matrix_signature_is_chunked_and_python_memory_bounded(monkeypatch):
    np = pytest.importorskip("numpy")
    sp = pytest.importorskip("scipy.sparse")
    n_rows = 100_000
    columns = np.array([0, 17, 1024, 2047], dtype=np.int32)
    indices = np.tile(columns, n_rows)
    indptr = np.arange(0, indices.size + 1, columns.size, dtype=np.int64)
    data = np.ones(indices.size, dtype=np.float32)
    matrix = sp.csr_matrix((data, indices, indptr), shape=(n_rows, 2048))

    def forbid_coo(*_args, **_kwargs):
        raise AssertionError("production-scale comparison must not materialise COO")

    monkeypatch.setattr(sp.csr_matrix, "tocoo", forbid_coo)
    tracemalloc.start()
    try:
        signature = compare_outputs._matrix_signature(matrix, require_finite=True)
        _current, peak = tracemalloc.get_traced_memory()
    finally:
        tracemalloc.stop()

    assert signature["nnz"] == 400_000
    assert "entries" not in signature
    assert peak < 32 * 1024 * 1024


@pytest.mark.integration
def test_h5ad_signature_rejects_integer_duplicate_sum_overflow():
    np = pytest.importorskip("numpy")
    sp = pytest.importorskip("scipy.sparse")
    matrix = sp.csr_matrix(
        (
            np.asarray([120, 120], dtype=np.int8),
            np.asarray([0, 0], dtype=np.int32),
            np.asarray([0, 2], dtype=np.int32),
        ),
        shape=(1, 1),
    )

    with pytest.raises(ValueError, match="integer overflow"):
        compare_outputs._matrix_signature(matrix, require_finite=True)


@pytest.mark.integration
def test_h5ad_accepts_valid_zero_observation_objects(tmp_path):
    anndata = pytest.importorskip("anndata")
    np = pytest.importorskip("numpy")
    pd = pytest.importorskip("pandas")
    sp = pytest.importorskip("scipy.sparse")
    path_a = tmp_path / "a.h5ad"
    path_b = tmp_path / "b.h5ad"
    for path in (path_a, path_b):
        adata = anndata.AnnData(
            sp.csr_matrix((0, 2), dtype=np.float32),
            obs=pd.DataFrame(index=pd.Index([], dtype=str)),
            var=pd.DataFrame(index=["Gene One", "Gene Two"]),
        )
        adata.write_h5ad(path)

    verdict, detail = compare_outputs._compare_h5ad(path_a, path_b, 10)

    assert verdict == compare_outputs.EQUAL
    assert detail == {"shape": [0, 2], "nnz": 0}


@pytest.mark.integration
@pytest.mark.parametrize(
    "invalid_kind", ["empty", "arbitrary", "corrupt", "truncated"]
)
def test_h5ad_rejects_invalid_files(tmp_path, invalid_kind):
    valid = tmp_path / "valid.h5ad"
    invalid = tmp_path / "invalid.h5ad"
    _make_adata([[1, 0], [2, 3]]).write_h5ad(valid)
    if invalid_kind == "empty":
        invalid.write_bytes(b"")
    elif invalid_kind == "arbitrary":
        invalid.write_bytes(b"these are not HDF5 bytes")
    elif invalid_kind == "corrupt":
        data = bytearray(valid.read_bytes())
        data[:8] = b"BROKEN!!"
        invalid.write_bytes(data)
    else:
        invalid.write_bytes(valid.read_bytes()[:-64])

    verdict, detail = compare_outputs._compare_h5ad(valid, invalid, 10)

    assert verdict == compare_outputs.DIFFER
    assert "b" in detail["validation_errors"]


@pytest.mark.unit
def test_h5ad_accepts_only_symmetric_named_zero_byte_sentinels(tmp_path):
    path_a, path_b = _paired_paths(tmp_path, "sample.raw_feature_bc_matrix.empty.h5ad")
    path_a.write_bytes(b"")
    path_b.write_bytes(b"")

    verdict, detail = compare_outputs._compare_h5ad(path_a, path_b, 10)

    assert verdict == compare_outputs.EQUAL
    assert detail["sentinel"] == "named empty H5AD"

    path_b.write_bytes(b"not empty")
    verdict, detail = compare_outputs._compare_h5ad(path_a, path_b, 10)
    assert verdict == compare_outputs.DIFFER
    assert "b" in detail["validation_errors"]


@pytest.mark.unit
def test_h5ad_rejects_symmetric_generic_zero_byte_files(tmp_path):
    path_a, path_b = _paired_paths(tmp_path, "sample.h5ad")
    path_a.write_bytes(b"")
    path_b.write_bytes(b"")

    verdict, detail = compare_outputs._compare_h5ad(path_a, path_b, 10)

    assert verdict == compare_outputs.DIFFER
    assert set(detail["validation_errors"]) == {"a", "b"}


@pytest.mark.unit
def test_h5ad_fails_closed_without_anndata(tmp_path, monkeypatch):
    path_a = tmp_path / "a.h5ad"
    path_b = tmp_path / "b.h5ad"
    path_a.write_bytes(b"same bytes")
    path_b.write_bytes(b"same bytes")
    monkeypatch.setattr(compare_outputs, "_HAVE_ANNDATA", False)

    verdict, detail = compare_outputs._compare_h5ad(path_a, path_b, 10)

    assert verdict == compare_outputs.DIFFER
    assert "dependency_error" in detail


SAMTOOLS = shutil.which("samtools")


def _write_bam(path, records=(), extra_header=(), sort_order="coordinate"):
    if SAMTOOLS is None:
        pytest.skip("samtools is not installed")
    sam = path.with_suffix(".sam")
    hd = "@HD\tVN:1.6"
    if sort_order is not None:
        hd += f"\tSO:{sort_order}"
    lines = [hd, "@SQ\tSN:chr1\tLN:1000"]
    lines.extend(extra_header)
    lines.extend(records)
    sam.write_text("\n".join(lines) + "\n", encoding="utf-8")
    subprocess.run(
        [SAMTOOLS, "view", "--no-PG", "-b", "-o", str(path), str(sam)],
        check=True,
        capture_output=True,
        text=True,
    )
    return path


def _sam_record(
    name="read1",
    pos=10,
    seq="ACGT",
    qual="IIII",
    gene="gene-1",
    cell="AAAC",
    umi="TTTT",
    flag=0,
    mapq=60,
    cigar="4M",
    mate="*\t0\t0",
    extra_tags=(),
    identity_tags=True,
    include_xt=True,
):
    tags = list(extra_tags)
    if identity_tags:
        tags.append(f"CB:Z:{cell}")
        if umi is not None:
            tags.append(f"UB:Z:{umi}")
    if include_xt:
        tags.append(f"XT:Z:{gene}")
    return (
        f"{name}\t{flag}\tchr1\t{pos}\t{mapq}\t{cigar}\t{mate}\t{seq}\t{qual}\t"
        + "\t".join(tags)
    )


@pytest.mark.unit
@pytest.mark.parametrize(
    ("relative", "stage", "contract", "require_xt"),
    [
        ("STAR/sample_Aligned.out.bam", "star", "alignment", False),
        (
            "STAR/sample_Aligned.sortedByCoord.out.bam",
            "star-empty",
            "alignment",
            False,
        ),
        (
            "featureCounts/sample_Aligned.sortedByCoord.out.bam.featureCounts.bam",
            "featurecounts",
            "alignment",
            False,
        ),
        (
            "featureCounts/sample.mapped.sorted.filtered.annotated.bam",
            "annotated",
            "alignment",
            True,
        ),
        (
            "deduplication/sample.dedup.bam",
            "dedup",
            "dedup-representative",
            True,
        ),
        ("custom/sample.bam", "generic", "alignment", False),
    ],
)
def test_bam_stage_contract_is_narrow_and_path_specific(
    tmp_path, relative, stage, contract, require_xt
):
    path = tmp_path / relative
    assert compare_outputs._bam_stage(path) == stage
    assert compare_outputs._bam_contract(path) == contract
    assert compare_outputs._bam_requires_xt(path) is require_xt


@pytest.mark.integration
def test_dedup_bam_normalises_only_representative_read_volatility(tmp_path):
    dir_a = tmp_path / "a" / "deduplication"
    dir_b = tmp_path / "b" / "deduplication"
    dir_a.mkdir(parents=True)
    dir_b.mkdir(parents=True)
    path_a = _write_bam(dir_a / "sample.dedup.bam", [_sam_record()])
    path_b = _write_bam(
        dir_b / "sample.dedup.bam",
        [_sam_record(
            name="another-source-read",
            seq="TGCA",
            qual="####",
            mapq=1,
            cigar="2M2S",
            mate="=\t99\t20",
            extra_tags=("AS:i:2", "NM:i:1"),
        )],
        extra_header=["@PG\tID:volatile\tPN:tool\tCL:tool --threads 2"],
    )

    verdict, detail = compare_outputs._compare_bam(path_a, path_b, 10)

    assert verdict == compare_outputs.EQUAL
    assert detail == {"count": 1}


@pytest.mark.integration
@pytest.mark.parametrize(
    ("base", "alternate"),
    [
        (
            {"pos": 10, "cigar": "4M"},
            {"pos": 12, "cigar": "2S2M"},
        ),
        (
            {"flag": 16, "pos": 10, "cigar": "4M"},
            {"flag": 16, "pos": 10, "cigar": "2M2S"},
        ),
    ],
)
def test_dedup_bam_uses_umi_tools_adjusted_five_prime_position(
    tmp_path, base, alternate
):
    dir_a = tmp_path / "a" / "deduplication"
    dir_b = tmp_path / "b" / "deduplication"
    dir_a.mkdir(parents=True)
    dir_b.mkdir(parents=True)
    path_a = _write_bam(dir_a / "sample.dedup.bam", [_sam_record(**base)])
    path_b = _write_bam(
        dir_b / "sample.dedup.bam", [_sam_record(**alternate)]
    )

    verdict, detail = compare_outputs._compare_bam(path_a, path_b, 10)

    assert verdict == compare_outputs.EQUAL
    assert detail == {"count": 1}


@pytest.mark.integration
def test_dedup_bam_external_sorts_adjusted_positions(tmp_path):
    """Adjusted 5' keys need not follow the BAM's raw coordinate order."""
    records_a = [
        _sam_record(name="plain", pos=7, cell="TTTG"),
        _sam_record(
            name="soft-clipped",
            pos=10,
            cigar="4S4M",
            seq="AAAACCCC",
            qual="IIIIIIII",
            cell="AAAC",
        ),
    ]
    records_b = [
        _sam_record(name="other-soft-source", pos=6, cell="AAAC"),
        _sam_record(
            name="other-plain-source",
            pos=9,
            cigar="2S4M",
            seq="AACCCC",
            qual="IIIIII",
            cell="TTTG",
        ),
    ]
    dir_a = tmp_path / "a" / "deduplication"
    dir_b = tmp_path / "b" / "deduplication"
    dir_a.mkdir(parents=True)
    dir_b.mkdir(parents=True)
    path_a = _write_bam(dir_a / "sample.dedup.bam", records_a)
    path_b = _write_bam(dir_b / "sample.dedup.bam", records_b)

    verdict, detail = compare_outputs._compare_bam(path_a, path_b, 10)

    assert verdict == compare_outputs.EQUAL
    assert detail == {"count": 2}


@pytest.mark.integration
def test_dedup_bam_rejects_changed_adjusted_five_prime_position(tmp_path):
    dir_a = tmp_path / "a" / "deduplication"
    dir_b = tmp_path / "b" / "deduplication"
    dir_a.mkdir(parents=True)
    dir_b.mkdir(parents=True)
    path_a = _write_bam(
        dir_a / "sample.dedup.bam", [_sam_record(pos=10, cigar="4M")]
    )
    path_b = _write_bam(
        dir_b / "sample.dedup.bam", [_sam_record(pos=11, cigar="4M")]
    )

    verdict, detail = compare_outputs._compare_bam(path_a, path_b, 10)

    assert verdict == compare_outputs.DIFFER
    assert detail["record_diffs"]


@pytest.mark.integration
def test_dedup_bam_accepts_read_name_identity_and_coordinate_tie_reordering(tmp_path):
    records_a = [
        _sam_record(name="source-a_AAAC_TTTT", identity_tags=False),
        _sam_record(name="source-b", cell="TTTG", umi=None),
    ]
    records_b = [
        _sam_record(name="different-source", cell="TTTG", umi=None),
        _sam_record(name="another-source_AAAC_TTTT", identity_tags=False),
    ]
    dir_a = tmp_path / "a" / "deduplication"
    dir_b = tmp_path / "b" / "deduplication"
    dir_a.mkdir(parents=True)
    dir_b.mkdir(parents=True)
    path_a = _write_bam(dir_a / "sample.dedup.bam", records_a)
    path_b = _write_bam(dir_b / "sample.dedup.bam", records_b)

    verdict, detail = compare_outputs._compare_bam(path_a, path_b, 10)

    assert verdict == compare_outputs.EQUAL
    assert detail == {"count": 2}


@pytest.mark.integration
def test_bam_normalises_published_unsorted_star_records_with_absent_xt(tmp_path):
    records_a = [
        _sam_record(name="source-a_AAAC_", pos=20, identity_tags=False, include_xt=False),
        _sam_record(name="source-b_TTTG_", pos=10, identity_tags=False, include_xt=False),
    ]
    records_b = list(reversed(records_a))
    dir_a = tmp_path / "a" / "STAR"
    dir_b = tmp_path / "b" / "STAR"
    dir_a.mkdir(parents=True)
    dir_b.mkdir(parents=True)
    path_a = _write_bam(
        dir_a / "sample_Aligned.out.bam", records_a, sort_order=None
    )
    path_b = _write_bam(
        dir_b / "sample_Aligned.out.bam", records_b, sort_order="unsorted"
    )

    verdict, detail = compare_outputs._compare_bam(path_a, path_b, 10)

    assert verdict == compare_outputs.EQUAL
    assert detail == {"count": 2}


@pytest.mark.integration
def test_only_configured_star_output_may_omit_sort_order(tmp_path):
    dir_a = tmp_path / "a" / "featureCounts"
    dir_b = tmp_path / "b" / "featureCounts"
    dir_a.mkdir(parents=True)
    dir_b.mkdir(parents=True)
    name = "sample_Aligned.sortedByCoord.out.bam.featureCounts.bam"
    record = _sam_record(include_xt=False)
    path_a = _write_bam(dir_a / name, [record], sort_order=None)
    path_b = _write_bam(dir_b / name, [record], sort_order=None)

    verdict, detail = compare_outputs._compare_bam(path_a, path_b, 10)

    assert verdict == compare_outputs.DIFFER
    assert set(detail["validation_errors"]) == {"a", "b"}
    assert "exactly one SO tag" in detail["validation_errors"]["a"]


@pytest.mark.integration
def test_alignment_comparison_does_not_trust_stale_coordinate_header(tmp_path):
    records = [
        _sam_record(name="later", pos=20),
        _sam_record(name="earlier", pos=10),
    ]
    dir_a = tmp_path / "a" / "featureCounts"
    dir_b = tmp_path / "b" / "featureCounts"
    dir_a.mkdir(parents=True)
    dir_b.mkdir(parents=True)
    filename = "sample.mapped.sorted.filtered.annotated.bam"
    path_a = _write_bam(dir_a / filename, records, sort_order="coordinate")
    path_b = _write_bam(
        dir_b / filename, list(reversed(records)), sort_order="coordinate"
    )

    verdict, detail = compare_outputs._compare_bam(path_a, path_b, 10)

    assert verdict == compare_outputs.EQUAL
    assert detail == {"count": 2}


@pytest.mark.integration
def test_star_empty_fallback_must_have_zero_records_on_both_sides(tmp_path):
    dir_a = tmp_path / "a" / "STAR"
    dir_b = tmp_path / "b" / "STAR"
    dir_a.mkdir(parents=True)
    dir_b.mkdir(parents=True)
    filename = "sample_Aligned.sortedByCoord.out.bam"
    record = _sam_record(include_xt=False)
    path_a = _write_bam(dir_a / filename, [record])
    path_b = _write_bam(dir_b / filename, [record])

    verdict, detail = compare_outputs._compare_bam(path_a, path_b, 10)

    assert verdict == compare_outputs.DIFFER
    assert set(detail["validation_errors"]) == {"a", "b"}
    assert "must contain zero" in detail["validation_errors"]["a"]


@pytest.mark.integration
@pytest.mark.parametrize(
    ("directory", "filename", "include_xt"),
    [
        ("STAR", "sample_Aligned.out.bam", False),
        (
            "featureCounts",
            "sample_Aligned.sortedByCoord.out.bam.featureCounts.bam",
            False,
        ),
        (
            "featureCounts",
            "sample.mapped.sorted.filtered.annotated.bam",
            True,
        ),
    ],
)
@pytest.mark.parametrize(
    ("difference", "changes"),
    [
        ("QNAME", {"name": "other-read"}),
        ("FLAG", {"flag": 256}),
        ("MAPQ", {"mapq": 1}),
        ("CIGAR", {"cigar": "2M2S"}),
        ("mate fields", {"mate": "=\t10\t4"}),
        ("SEQ", {"seq": "TGCA"}),
        ("QUAL", {"qual": "####"}),
        ("optional tags", {"extra_tags": ("AS:i:2",)}),
    ],
)
def test_pre_dedup_bams_reject_changed_alignment_content(
    tmp_path, directory, filename, include_xt, difference, changes
):
    dir_a = tmp_path / "a" / directory
    dir_b = tmp_path / "b" / directory
    dir_a.mkdir(parents=True)
    dir_b.mkdir(parents=True)
    base_record = _sam_record(include_xt=include_xt)
    changed_record = _sam_record(include_xt=include_xt, **changes)
    sort_order = "unsorted" if directory == "STAR" else "coordinate"
    path_a = _write_bam(dir_a / filename, [base_record], sort_order=sort_order)
    path_b = _write_bam(dir_b / filename, [changed_record], sort_order=sort_order)

    verdict, detail = compare_outputs._compare_bam(path_a, path_b, 10)

    assert verdict == compare_outputs.DIFFER, difference
    assert detail["record_diffs"], difference


@pytest.mark.integration
def test_pre_dedup_bam_ignores_only_optional_tag_order_and_pg_header(tmp_path):
    dir_a = tmp_path / "a" / "featureCounts"
    dir_b = tmp_path / "b" / "featureCounts"
    dir_a.mkdir(parents=True)
    dir_b.mkdir(parents=True)
    name = "sample_Aligned.sortedByCoord.out.bam.featureCounts.bam"
    tags_a = ("AS:i:2", "NM:i:1")
    tags_b = tuple(reversed(tags_a))
    path_a = _write_bam(
        dir_a / name,
        [_sam_record(extra_tags=tags_a)],
        extra_header=["@PG\tID:run-a\tPN:tool\tCL:tool --threads 1"],
    )
    path_b = _write_bam(
        dir_b / name,
        [_sam_record(extra_tags=tags_b)],
        extra_header=["@PG\tID:run-b\tPN:tool\tCL:tool --threads 8"],
    )

    verdict, detail = compare_outputs._compare_bam(path_a, path_b, 10)

    assert verdict == compare_outputs.EQUAL
    assert detail == {"count": 1}


@pytest.mark.integration
@pytest.mark.parametrize(
    ("directory", "filename"),
    [
        (
            "featureCounts",
            "sample.mapped.sorted.filtered.annotated.bam",
        ),
        ("deduplication", "sample.dedup.bam"),
    ],
)
def test_bam_xt_presence_contract_follows_published_bam_stage(
    tmp_path, directory, filename
):
    # The initial featureCounts BAM contains legitimate unassigned records with
    # no XT, so absence remains an explicit stable empty key at that stage.
    initial_a = tmp_path / "initial-a" / "featureCounts"
    initial_b = tmp_path / "initial-b" / "featureCounts"
    initial_a.mkdir(parents=True)
    initial_b.mkdir(parents=True)
    name = "sample_Aligned.sortedByCoord.out.bam.featureCounts.bam"
    record = _sam_record(include_xt=False)
    path_a = _write_bam(initial_a / name, [record])
    path_b = _write_bam(initial_b / name, [record])
    assert compare_outputs._compare_bam(path_a, path_b, 10)[0] \
        == compare_outputs.EQUAL

    # High-confidence and deduplicated BAMs contain only gene-assigned records;
    # symmetric XT loss there is invalid rather than equivalent.
    required_a = tmp_path / "required-a" / directory
    required_b = tmp_path / "required-b" / directory
    required_a.mkdir(parents=True)
    required_b.mkdir(parents=True)
    path_a = _write_bam(required_a / filename, [record])
    path_b = _write_bam(required_b / filename, [record])
    verdict, detail = compare_outputs._compare_bam(path_a, path_b, 10)
    assert verdict == compare_outputs.DIFFER
    assert set(detail["validation_errors"]) == {"a", "b"}
    assert "must contain an XT" in detail["validation_errors"]["a"]


@pytest.mark.integration
def test_bam_stable_molecule_multiplicity_is_part_of_contract(tmp_path):
    dir_a = tmp_path / "a" / "deduplication"
    dir_b = tmp_path / "b" / "deduplication"
    dir_a.mkdir(parents=True)
    dir_b.mkdir(parents=True)
    path_a = _write_bam(
        dir_a / "sample.dedup.bam",
        [
            _sam_record(name="a1", cell="AAAC"),
            _sam_record(name="a2", cell="AAAC"),
            _sam_record(name="a3", cell="TTTG"),
        ],
    )
    path_b = _write_bam(
        dir_b / "sample.dedup.bam",
        [
            _sam_record(name="b1", cell="AAAC"),
            _sam_record(name="b2", cell="TTTG"),
            _sam_record(name="b3", cell="TTTG"),
        ],
    )

    verdict, detail = compare_outputs._compare_bam(path_a, path_b, 10)

    assert verdict == compare_outputs.DIFFER
    assert detail["record_diffs"][0]["only_a"][0]["count"] == 1
    assert detail["record_diffs"][0]["only_b"][0]["count"] == 1


@pytest.mark.integration
@pytest.mark.parametrize("difference", ["position", "gene", "cell", "umi", "strand"])
def test_bam_rejects_same_count_and_flagstat_with_different_alignments(
    tmp_path, difference
):
    dir_a = tmp_path / "a" / "deduplication"
    dir_b = tmp_path / "b" / "deduplication"
    dir_a.mkdir(parents=True)
    dir_b.mkdir(parents=True)
    path_a = _write_bam(dir_a / "sample.dedup.bam", [_sam_record()])
    changes = {
        "position": {"pos": 20},
        "gene": {"gene": "gene-2"},
        "cell": {"cell": "TTTG"},
        "umi": {"umi": "GGGG"},
        "strand": {"flag": 16},
    }
    record_b = _sam_record(**changes[difference])
    path_b = _write_bam(dir_b / "sample.dedup.bam", [record_b])

    verdict, detail = compare_outputs._compare_bam(path_a, path_b, 10)

    assert verdict == compare_outputs.DIFFER
    assert detail["record_diffs"]


@pytest.mark.integration
def test_unknown_bam_path_uses_conservative_complete_alignment_contract(tmp_path):
    dir_a = tmp_path / "a" / "custom"
    dir_b = tmp_path / "b" / "custom"
    dir_a.mkdir(parents=True)
    dir_b.mkdir(parents=True)
    path_a = _write_bam(dir_a / "sample.bam", [_sam_record()])
    path_b = _write_bam(
        dir_b / "sample.bam", [_sam_record(cigar="2M2S", seq="TGCA")]
    )

    verdict, detail = compare_outputs._compare_bam(path_a, path_b, 10)

    assert verdict == compare_outputs.DIFFER
    assert detail["record_diffs"]


@pytest.mark.integration
def test_bam_accepts_valid_header_only_files(tmp_path):
    path_a = _write_bam(tmp_path / "a.bam")
    path_b = _write_bam(tmp_path / "b.bam")

    verdict, detail = compare_outputs._compare_bam(path_a, path_b, 10)

    assert verdict == compare_outputs.EQUAL
    assert detail == {"count": 0}


@pytest.mark.integration
@pytest.mark.parametrize(
    "invalid_kind", ["empty", "arbitrary", "corrupt", "truncated"]
)
def test_bam_rejects_invalid_files(tmp_path, invalid_kind):
    valid = _write_bam(tmp_path / "valid.bam", [_sam_record()])
    invalid = tmp_path / "invalid.bam"
    if invalid_kind == "empty":
        invalid.write_bytes(b"")
    elif invalid_kind == "arbitrary":
        invalid.write_bytes(b"not a BAM")
    elif invalid_kind == "corrupt":
        data = bytearray(valid.read_bytes())
        first_block_size = int.from_bytes(data[16:18], "little") + 1
        compressed_payload_start = 18
        compressed_payload_end = first_block_size - 8
        corrupt_at = (compressed_payload_start + compressed_payload_end) // 2
        data[corrupt_at] ^= 0xFF
        invalid.write_bytes(data)
    else:
        invalid.write_bytes(valid.read_bytes()[:-16])

    verdict, detail = compare_outputs._compare_bam(valid, invalid, 10)

    assert verdict == compare_outputs.DIFFER
    assert "b" in detail["validation_errors"]


@pytest.mark.unit
def test_bam_fails_closed_without_samtools(tmp_path, monkeypatch):
    path_a = tmp_path / "a.bam"
    path_b = tmp_path / "b.bam"
    path_a.write_bytes(b"same bytes")
    path_b.write_bytes(b"same bytes")
    monkeypatch.setattr(compare_outputs, "_HAVE_SAMTOOLS", False)

    verdict, detail = compare_outputs._compare_bam(path_a, path_b, 10)

    assert verdict == compare_outputs.DIFFER
    assert "dependency_error" in detail
