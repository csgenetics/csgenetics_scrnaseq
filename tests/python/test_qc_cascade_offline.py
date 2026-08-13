"""Offline and workflow-contract tests for the published QC cascade plots.

The consolidated report consumes small Plotly fragments, but the HTML files
published alongside the CSV metrics are opened directly by customers. Those
published files therefore need their own inline Plotly runtime. These tests
keep the two output classes separate and render both standalone variants in an
offline browser when Playwright/Chromium are available.
"""

from html.parser import HTMLParser
import csv
import math
import os
from pathlib import Path
import subprocess
import sys

import pytest


pytestmark = pytest.mark.integration

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT = REPO_ROOT / "bin" / "qc_cascade_plot.py"
FIXTURE_DIR = Path(__file__).parent / "fixtures" / "report"
sys.path.insert(0, os.fspath(REPO_ROOT / "bin"))
import create_consolidated_report as consolidated_report  # noqa: E402
import qc_cascade_plot as qc_cascade  # noqa: E402


class _ExternalResourceParser(HTMLParser):
    """Collect active remote script/link/image resources, not inert JS text."""

    def __init__(self):
        super().__init__()
        self.remote_resources = []

    def handle_starttag(self, tag, attrs):
        attribute = "href" if tag == "link" else "src"
        if tag not in {"link", "script", "img", "iframe"}:
            return
        value = dict(attrs).get(attribute, "")
        if value.startswith(("http://", "https://", "//")):
            self.remote_resources.append(value)


@pytest.fixture(scope="module")
def qc_outputs(tmp_path_factory):
    work = tmp_path_factory.mktemp("qc-cascade-offline")
    sample_1 = FIXTURE_DIR / "SAMPLE1.metrics.csv"
    sample_2 = FIXTURE_DIR / "SAMPLE2.metrics.csv"

    single = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--mode",
            "single",
            "--sample-id",
            "SAMPLE1",
            "--metrics-csv",
            str(sample_1),
        ],
        cwd=work,
        capture_output=True,
        text=True,
    )
    assert single.returncode == 0, single.stderr

    multi = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--mode",
            "multi",
            "--csv-files",
            str(sample_1),
            str(sample_2),
        ],
        cwd=work,
        capture_output=True,
        text=True,
    )
    assert multi.returncode == 0, multi.stderr

    outputs = {
        "single_fragment": work / "SAMPLE1.qc_cascade.html",
        "single_standalone": work / "SAMPLE1.qc_cascade.standalone.html",
        "multi_fragment": work / "multisample_qc_cascade.fragment.html",
        "multi_standalone": work / "multisample_qc_cascade.html",
    }
    for name, output in outputs.items():
        assert output.is_file(), f"missing {name}: {output}"
    return outputs


def _read(path):
    return path.read_text(encoding="utf-8")


def _write_zero_metrics(path, sample_fixture="SAMPLE1.metrics.csv"):
    with (FIXTURE_DIR / sample_fixture).open(newline="", encoding="utf-8") as source:
        rows = list(csv.reader(source))
    for row in rows[1:]:
        row[1] = "0"
    with path.open("w", newline="", encoding="utf-8") as destination:
        csv.writer(destination, lineterminator="\n").writerows(rows)
    return path


def _all_figure_y_values_are_finite(fig):
    return all(
        math.isfinite(float(value))
        for trace in fig.data
        for value in (trace.y if trace.y is not None else ())
    )


def test_required_read_qc_schema_accepts_legitimate_zero_values(tmp_path):
    metrics_path = _write_zero_metrics(tmp_path / "EMPTY.metrics.csv")
    metrics = qc_cascade.QCCascadePlotter("single").read_metrics_csv(metrics_path)

    assert set(metrics) == set(qc_cascade.REQUIRED_READ_QC_METRICS)
    assert set(metrics.values()) == {0}


@pytest.mark.parametrize(
    "invalid_value", ["not-a-number", "NaN", "Infinity", "-1", "1.5"]
)
def test_required_metric_malformed_nonfinite_or_negative_fails(tmp_path, invalid_value):
    path = _write_zero_metrics(tmp_path / "EMPTY.metrics.csv")
    rows = list(csv.reader(path.open(newline="", encoding="utf-8")))
    rows[1][1] = invalid_value
    with path.open("w", newline="", encoding="utf-8") as destination:
        csv.writer(destination, lineterminator="\n").writerows(rows)

    with pytest.raises(ValueError, match="Required metric 'reads_pre_qc'"):
        qc_cascade.QCCascadePlotter("single").read_metrics_csv(path)


def test_integer_valued_float_count_is_accepted_and_normalized(tmp_path):
    path = _write_zero_metrics(tmp_path / "EMPTY.metrics.csv")
    rows = list(csv.reader(path.open(newline="", encoding="utf-8")))
    rows[1][1] = "1.0"
    with path.open("w", newline="", encoding="utf-8") as destination:
        csv.writer(destination, lineterminator="\n").writerows(rows)

    metrics = qc_cascade.QCCascadePlotter("single").read_metrics_csv(path)
    assert metrics["reads_pre_qc"] == 1
    assert isinstance(metrics["reads_pre_qc"], int)


def test_zero_input_cannot_mask_nonzero_downstream_count(tmp_path):
    path = _write_zero_metrics(tmp_path / "EMPTY.metrics.csv")
    rows = list(csv.reader(path.open(newline="", encoding="utf-8")))
    for row in rows:
        if row and row[0] == "polya_retained":
            row[1] = "1"
            break
    with path.open("w", newline="", encoding="utf-8") as destination:
        csv.writer(destination, lineterminator="\n").writerows(rows)

    with pytest.raises(
        ValueError, match="zero reads_pre_qc.*nonzero.*polya_retained"
    ):
        qc_cascade.QCCascadePlotter("single").read_metrics_csv(path)


def test_missing_required_metric_fails(tmp_path):
    path = _write_zero_metrics(tmp_path / "EMPTY.metrics.csv")
    rows = list(csv.reader(path.open(newline="", encoding="utf-8")))
    rows = [row for row in rows if not row or row[0] != "polya_retained"]
    with path.open("w", newline="", encoding="utf-8") as destination:
        csv.writer(destination, lineterminator="\n").writerows(rows)

    with pytest.raises(ValueError, match="missing required Read-QC metric.*polya_retained"):
        qc_cascade.QCCascadePlotter("single").read_metrics_csv(path)


def test_required_metric_with_wrong_classification_fails(tmp_path):
    path = _write_zero_metrics(tmp_path / "EMPTY.metrics.csv")
    rows = list(csv.reader(path.open(newline="", encoding="utf-8")))
    rows[1][4] = "Alignment QC"
    with path.open("w", newline="", encoding="utf-8") as destination:
        csv.writer(destination, lineterminator="\n").writerows(rows)

    with pytest.raises(ValueError, match="expected 'Read QC'"):
        qc_cascade.QCCascadePlotter("single").read_metrics_csv(path)


def test_invalid_metrics_header_fails(tmp_path):
    path = _write_zero_metrics(tmp_path / "EMPTY.metrics.csv")
    rows = list(csv.reader(path.open(newline="", encoding="utf-8")))
    rows[0][0] = "wrong_name"
    with path.open("w", newline="", encoding="utf-8") as destination:
        csv.writer(destination, lineterminator="\n").writerows(rows)

    with pytest.raises(ValueError, match="invalid header"):
        qc_cascade.QCCascadePlotter("single").read_metrics_csv(path)


def test_duplicate_required_metric_fails(tmp_path):
    path = _write_zero_metrics(tmp_path / "EMPTY.metrics.csv")
    rows = list(csv.reader(path.open(newline="", encoding="utf-8")))
    rows.append(rows[1])
    with path.open("w", newline="", encoding="utf-8") as destination:
        csv.writer(destination, lineterminator="\n").writerows(rows)

    with pytest.raises(ValueError, match="occurs more than once"):
        qc_cascade.QCCascadePlotter("single").read_metrics_csv(path)


def test_malformed_metrics_csv_row_fails(tmp_path):
    path = _write_zero_metrics(tmp_path / "EMPTY.metrics.csv")
    with path.open("a", encoding="utf-8") as destination:
        destination.write("broken,row\n")

    with pytest.raises(ValueError, match="has 2 columns; expected 5"):
        qc_cascade.QCCascadePlotter("single").read_metrics_csv(path)


def test_zero_input_single_and_multi_figures_have_only_finite_y_values(
    tmp_path, monkeypatch
):
    sample_1 = _write_zero_metrics(tmp_path / "EMPTY1.metrics.csv")
    sample_2 = _write_zero_metrics(
        tmp_path / "EMPTY2.metrics.csv", sample_fixture="SAMPLE2.metrics.csv"
    )
    monkeypatch.chdir(tmp_path)
    plotter = qc_cascade.QCCascadePlotter("single", sample_id="EMPTY1")
    single_figure = plotter.create_single_sample_plot(
        plotter.read_metrics_csv(sample_1)
    )
    multi_figure = qc_cascade.QCCascadePlotter("multi").create_multi_sample_plot(
        [sample_1, sample_2]
    )

    assert _all_figure_y_values_are_finite(single_figure)
    assert _all_figure_y_values_are_finite(multi_figure)
    assert "NaN" not in _read(tmp_path / "EMPTY1.qc_cascade.html")
    assert "Infinity" not in _read(tmp_path / "EMPTY1.qc_cascade.html")
    assert "NaN" not in _read(tmp_path / "multisample_qc_cascade.fragment.html")
    assert "Infinity" not in _read(tmp_path / "multisample_qc_cascade.fragment.html")


def test_two_empty_samples_generate_complete_consolidated_report(tmp_path):
    sample_ids = ("EMPTY1", "EMPTY2")
    metrics_paths = []
    for index, sample_id in enumerate(sample_ids, start=1):
        metrics = _write_zero_metrics(
            tmp_path / f"{sample_id}.metrics.csv",
            sample_fixture=f"SAMPLE{index}.metrics.csv",
        )
        metrics_paths.append(metrics)
        (tmp_path / f"{sample_id}_counts_pdf_with_threshold.html").touch()

        single = subprocess.run(
            [
                sys.executable,
                str(SCRIPT),
                "--mode",
                "single",
                "--sample-id",
                sample_id,
                "--metrics-csv",
                str(metrics),
            ],
            cwd=tmp_path,
            capture_output=True,
            text=True,
        )
        assert single.returncode == 0, single.stderr

    multi = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--mode",
            "multi",
            "--csv-files",
            *(str(path) for path in metrics_paths),
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
    )
    assert multi.returncode == 0, multi.stderr

    report = subprocess.run(
        [
            sys.executable,
            str(REPO_ROOT / "bin" / "create_consolidated_report.py"),
            str(REPO_ROOT / "templates" / "consolidated_report_template.html.jinja2"),
            "FALSE",
            str(REPO_ROOT / "assets" / "vendor"),
            str(tmp_path / "multisample_qc_cascade.fragment.html"),
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
    )

    assert report.returncode == 0, report.stderr
    html = _read(tmp_path / "consolidated_report.html")
    assert all(sample_id in html for sample_id in sample_ids)
    multisample_csv = tmp_path / "multisample_out.csv"
    assert multisample_csv.is_file()
    with multisample_csv.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.reader(handle))
    assert rows[0][-2:] == list(sample_ids)
    assert rows[1:]
    for row in rows[1:]:
        assert len(row) == 6
        values = [float(value) for value in row[-2:]]
        assert all(math.isfinite(value) and value == 0 for value in values)
    for fragment in (
        tmp_path / "EMPTY1.qc_cascade.html",
        tmp_path / "EMPTY2.qc_cascade.html",
        tmp_path / "multisample_qc_cascade.fragment.html",
    ):
        consolidated_report.validate_plot_fragment(_read(fragment))


@pytest.mark.parametrize("variant", ["single", "multi"])
def test_fragment_and_standalone_outputs_are_distinct(qc_outputs, variant):
    fragment = _read(qc_outputs[f"{variant}_fragment"])
    standalone = _read(qc_outputs[f"{variant}_standalone"])

    assert "Plotly.newPlot" in fragment
    assert "<html" not in fragment.lower()
    assert len(fragment) < 1_000_000, "report fragment unexpectedly bundles Plotly.js"

    assert "<!doctype html>" in standalone.lower()
    assert "Plotly.newPlot" in standalone
    assert len(standalone) > 1_000_000, "standalone output does not embed Plotly.js"

    parser = _ExternalResourceParser()
    parser.feed(standalone)
    assert parser.remote_resources == [], (
        f"standalone {variant} output references remote resources: "
        f"{parser.remote_resources}"
    )


def test_multi_sample_labels_do_not_include_staging_paths(qc_outputs):
    """Absolute staged input paths must not leak into customer-visible labels."""
    html = _read(qc_outputs["multi_standalone"])
    assert "SAMPLE1" in html
    assert "SAMPLE2" in html
    assert os.fspath(qc_outputs["multi_standalone"].parent) not in html


@pytest.mark.parametrize("output_name", ["single_fragment", "multi_fragment"])
def test_generated_fragment_satisfies_consolidated_report_contract(qc_outputs, output_name):
    raw = _read(qc_outputs[output_name])
    trusted = consolidated_report.validate_plot_fragment(raw)

    assert "Plotly.newPlot" in trusted
    assert f'nonce="{consolidated_report.REPORT_SCRIPT_NONCE}"' in trusted


def test_nextflow_keeps_internal_fragments_out_of_published_outputs():
    single_module = (REPO_ROOT / "modules/local/qc_cascade_plot_single/main.nf").read_text()
    multi_module = (REPO_ROOT / "modules/local/qc_cascade_plot_multi/main.nf").read_text()
    workflow = (REPO_ROOT / "main.nf").read_text()

    assert "emit: qc_cascade_fragment" in single_module
    assert "emit: qc_cascade_report" in single_module
    assert 'pattern: "*.qc_cascade.standalone.html"' in single_module
    assert 'saveAs: { "${sample_id}.qc_cascade.html" }' in single_module

    assert "emit: qc_cascade_fragment" in multi_module
    assert "emit: qc_cascade_report" in multi_module
    assert 'pattern: "multisample_qc_cascade.html"' in multi_module

    assert "qc_cascade_plot_single.out.qc_cascade_fragment" in workflow
    assert "qc_cascade_plot_multi.out.qc_cascade_fragment" in workflow
    assert "ch_all_qc_cascade_fragments" in workflow
    assert ".out.qc_cascade_report" not in workflow
    assert ".out.qc_cascade_plot" not in workflow


@pytest.fixture(scope="module")
def chromium_browser():
    try:
        from playwright.sync_api import sync_playwright
    except ImportError:
        pytest.skip("Playwright is not installed")

    playwright = None
    try:
        playwright = sync_playwright().start()
        browser = playwright.chromium.launch()
    except Exception as exc:
        if playwright is not None:
            playwright.stop()
        pytest.skip(f"Playwright Chromium is unavailable: {exc}")

    yield browser
    browser.close()
    playwright.stop()


@pytest.mark.parametrize("output_name", ["single_standalone", "multi_standalone"])
def test_standalone_plot_renders_offline_without_remote_requests_or_js_errors(
    qc_outputs, chromium_browser, output_name
):
    context = chromium_browser.new_context(offline=True)
    page = context.new_page()
    errors = []
    remote_requests = []
    page.on("console", lambda message: errors.append(message.text) if message.type == "error" else None)
    page.on("pageerror", lambda error: errors.append(str(error)))
    page.on(
        "request",
        lambda request: remote_requests.append(request.url)
        if request.url.startswith(("http://", "https://"))
        else None,
    )

    page.goto(qc_outputs[output_name].as_uri(), wait_until="load")
    page.wait_for_selector(".js-plotly-plot .main-svg", state="attached", timeout=10_000)

    assert page.evaluate("typeof window.Plotly !== 'undefined'")
    assert page.locator(".js-plotly-plot .main-svg").count() >= 1
    assert remote_requests == [], f"offline plot attempted remote requests: {remote_requests}"
    assert errors == [], f"JavaScript errors while rendering offline plot: {errors}"
    context.close()
