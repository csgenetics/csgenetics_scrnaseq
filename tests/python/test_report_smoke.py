"""
Headless-browser smoke test for the consolidated HTML report.

Two real report bugs shipped that file-level / unit checks could not catch -- both
only manifested when the report was actually RENDERED in a browser:

  * Plotly.js was inlined AFTER the inline ``Plotly.newPlot()`` fragment calls, so
    every plot threw "Plotly is not defined" and rendered blank (the metric tables,
    which are plain HTML, still rendered, masking the failure).
  * The sample dropdown's text input doubled as both the search box and the
    selected-value display, so after selecting a sample the displayed sample id was
    re-applied as a filter query and the option list collapsed to one entry.

This test generates a report from a small committed fixture, opens it in headless
Chromium via Playwright, and asserts the report is FUNCTIONAL, not merely present:
no JS console errors, ``window.Plotly`` defined, every plot actually drew an SVG,
the dropdown selects the right pane and keeps all samples on reopen, the print
stylesheet reveals every per-sample pane, and the re-emitted ``multisample_out.csv``
keeps raw (separator-free) numbers so it stays machine-readable.

It is intentionally dependency-light: the report generator needs only jinja2 +
plotly, so this runs in a small Playwright CI image without the heavy pipeline
Python environment.
"""

import os
import shutil
import subprocess
import sys

import pytest

# Skip cleanly where Playwright is not installed (e.g. a bare dev env); CI installs it.
pytest.importorskip("playwright.sync_api")
from playwright.sync_api import sync_playwright  # noqa: E402

# Runs the report generator as a subprocess and renders it in a browser.
pytestmark = pytest.mark.integration

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
GENERATOR = os.path.join(REPO_ROOT, "bin", "create_consolidated_report.py")
TEMPLATE = os.path.join(REPO_ROOT, "templates", "consolidated_report_template.html.jinja2")
VENDOR_DIR = os.path.join(REPO_ROOT, "assets", "vendor")
FIXTURE_DIR = os.path.join(os.path.dirname(__file__), "fixtures", "report")

EXPECTED_SAMPLES = ["SAMPLE1", "SAMPLE2"]
# fixture has, per sample: 1 cell-caller + 1 qc-cascade plot, plus 1 multi-sample
# qc-cascade across all samples => 2*2 + 1 = 5 Plotly figures.
EXPECTED_PLOTS = 5


@pytest.fixture(scope="module")
def report_html(tmp_path_factory):
    """Generate the consolidated report from the fixture; return the .html path."""
    work = tmp_path_factory.mktemp("report")
    for name in os.listdir(FIXTURE_DIR):
        shutil.copy(os.path.join(FIXTURE_DIR, name), work / name)

    # Same invocation Nextflow uses: <template> <mixed_species> <vendor_dir> <multi_qc_cascade>
    subprocess.run(
        [sys.executable, GENERATOR, TEMPLATE, "FALSE", VENDOR_DIR,
         str(work / "multisample_qc_cascade.html")],
        cwd=work, check=True,
    )
    html = work / "consolidated_report.html"
    assert html.exists(), "generator did not produce consolidated_report.html"
    return html


@pytest.fixture(scope="module")
def browser():
    with sync_playwright() as p:
        b = p.chromium.launch()
        yield b
        b.close()


@pytest.fixture
def page(browser, report_html):
    """A freshly-loaded page per test, with JS errors collected from load onward."""
    pg = browser.new_page()
    errors = []
    pg.on("console", lambda m: errors.append(m.text) if m.type == "error" else None)
    pg.on("pageerror", lambda e: errors.append(str(e)))
    pg.goto(report_html.as_uri())
    pg.wait_for_timeout(1200)  # let inline Plotly.newPlot() calls settle
    pg._console_errors = errors
    yield pg
    pg.close()


def test_no_js_console_errors(page):
    assert page._console_errors == [], f"JS errors in rendered report: {page._console_errors}"


def test_plotly_is_defined(page):
    assert page.evaluate("typeof window.Plotly !== 'undefined'"), \
        "window.Plotly is undefined (script-ordering / missing-library regression)"


def test_every_plot_drew_an_svg(page):
    n_divs = page.eval_on_selector_all(".js-plotly-plot", "els => els.length")
    n_with_svg = page.eval_on_selector_all(
        ".js-plotly-plot", "els => els.filter(e => e.querySelector('.main-svg')).length")
    assert n_divs >= EXPECTED_PLOTS, f"expected >= {EXPECTED_PLOTS} plot divs, got {n_divs}"
    assert n_with_svg == n_divs, \
        f"{n_divs - n_with_svg} of {n_divs} plots drew no SVG (blank-plot regression)"


def test_dropdown_selects_the_correct_pane(page):
    page.click("#samplePickerInput")
    page.click(".cs-sample-option >> nth=1")  # select the 2nd sample
    shown = page.eval_on_selector_all(".cs-sample-pane.cs-show", "els => els.map(e => e.id)")
    assert shown == ["sample-pane-2"], f"expected only sample-pane-2 shown, got {shown}"


def test_dropdown_keeps_all_samples_after_select_and_reopen(page):
    page.click("#samplePickerInput")
    page.click(".cs-sample-option >> nth=1")  # select
    page.click("#samplePickerInput")          # reopen
    visible = page.eval_on_selector_all(
        ".cs-sample-option", "els => els.filter(e => !e.classList.contains('cs-hidden')).length")
    assert visible == len(EXPECTED_SAMPLES), \
        f"dropdown collapsed to {visible} option(s) after select+reopen (collapse regression)"


def test_dropdown_is_keyboard_operable(page):
    """The sample picker is a role=combobox; arrow keys must move the keyboard
    cursor (aria-activedescendant) and Enter must select the highlighted option."""
    page.focus("#samplePickerInput")
    page.keyboard.press("ArrowDown")
    page.keyboard.press("ArrowDown")
    ad = page.eval_on_selector("#samplePickerInput", "e => e.getAttribute('aria-activedescendant')")
    assert ad == "sample-opt-2", f"ArrowDown did not move the keyboard cursor (aria-activedescendant={ad})"
    page.keyboard.press("Enter")
    shown = page.eval_on_selector_all(".cs-sample-pane.cs-show", "els => els.map(e => e.id)")
    assert shown == ["sample-pane-2"], f"Enter did not select the highlighted sample, got {shown}"


def test_accessibility_structure(page):
    """Landmarks, a single page heading, a skip link, and focusable (not hover-only)
    info affordances -- the structural a11y guarantees."""
    assert page.eval_on_selector_all("h1", "els => els.length") == 1, "expected exactly one <h1>"
    assert page.query_selector("main#main-content") is not None, "missing <main> landmark"
    assert page.query_selector("a.cs-skip-link") is not None, "missing skip-to-content link"
    assert page.eval_on_selector_all("i.cs-info", "els => els.length") == 0, \
        "metric-info affordances must be focusable <button>s, not hover-only <i>"
    assert page.eval_on_selector_all("button.cs-info", "els => els.length") > 0, \
        "expected focusable info buttons"


def test_print_media_reveals_all_panes(page):
    page.emulate_media(media="print")
    visible = page.eval_on_selector_all(
        ".cs-sample-pane", "els => els.filter(e => e.offsetParent !== null).length")
    page.emulate_media(media="screen")
    assert visible == len(EXPECTED_SAMPLES), \
        f"print should reveal all {len(EXPECTED_SAMPLES)} panes (else PDF drops samples), got {visible}"


def test_csv_keeps_raw_separatorless_numbers(report_html):
    """The on-screen tables get thousands separators, but multisample_out.csv must
    stay machine-readable -- commas would corrupt the comma-delimited data file."""
    csv_path = os.path.join(os.path.dirname(report_html), "multisample_out.csv")
    with open(csv_path) as fh:
        for line in fh:
            if line.startswith("reads_pre_qc,"):
                values = line.rstrip("\n").split(",")[4:]
                assert values, "reads_pre_qc row has no sample values"
                for v in values:
                    assert v.isdigit(), f"CSV integer value carries a separator or is non-numeric: {v!r}"
                return
    pytest.fail("reads_pre_qc row not found in multisample_out.csv")


def test_run_provenance_renders_when_supplied(tmp_path):
    """When the consolidated_report process passes a provenance JSON, the report
    shows a Run-provenance header, a Key-results table, an Output-files pointer and a
    run id in the title and footer. (The fixture-based tests above pass NO provenance,
    confirming the sections are cleanly hidden when absent.)"""
    import json
    for name in os.listdir(FIXTURE_DIR):
        shutil.copy(os.path.join(FIXTURE_DIR, name), tmp_path / name)
    prov = json.dumps({
        "genome": "GRCh38", "annotation": "gencode.v44.gtf", "mixed": False,
        "pipeline_ver": "2.0.0", "commit": "abc1234", "revision": "main",
        "run_name": "cheeky_curie", "session_id": "e495-3a25", "start": "2026-06-17T09:00:00Z",
        "nf_version": "26.04.1", "outdir": "s3://example-bucket/run42",
        "barcode_kit": "IDT_IO_kit_v2.csv", "count_threshold": 100,
        "homepage": "https://github.com/csgenetics/csgenetics_scrnaseq",
    })
    subprocess.run(
        [sys.executable, GENERATOR, TEMPLATE, "FALSE", VENDOR_DIR,
         str(tmp_path / "multisample_qc_cascade.html"), prov],
        cwd=tmp_path, check=True,
    )
    html = (tmp_path / "consolidated_report.html").read_text()
    for token in ["Run provenance", "GRCh38", "cheeky_curie", "Key results", "Output files"]:
        assert token in html, f"provenance render missing {token!r}"
    assert "<title>CS Genetics scRNA-seq report - cheeky_curie" in html, "title lacks run id"


def test_malicious_sample_id_is_escaped_not_executed(tmp_path, browser):
    """Sample ids come from the customer's input sheet -- the one piece of
    untrusted text in the report. The generator runs with Jinja autoescaping on,
    so an HTML/JS payload in a sample id must be rendered as inert text, not
    injected as markup or executed as script."""
    payload = "<img src=x onerror=window.__xss=1>"
    for name in os.listdir(FIXTURE_DIR):
        dst = name.replace("SAMPLE1", payload) if name.startswith("SAMPLE1") else name
        shutil.copy(os.path.join(FIXTURE_DIR, name), tmp_path / dst)
    subprocess.run(
        [sys.executable, GENERATOR, TEMPLATE, "FALSE", VENDOR_DIR,
         str(tmp_path / "multisample_qc_cascade.html")],
        cwd=tmp_path, check=True,
    )
    html = (tmp_path / "consolidated_report.html").read_text()
    assert payload not in html, "raw (unescaped) sample-id payload was injected into the HTML"
    assert "&lt;img src=x onerror=window.__xss=1&gt;" in html, "payload was not present in escaped form"

    # And confirm it does not execute when rendered.
    pg = browser.new_page()
    pg.goto((tmp_path / "consolidated_report.html").as_uri())
    pg.wait_for_timeout(800)
    assert pg.evaluate("window.__xss === undefined"), "sample-id payload executed as script (XSS)"
    pg.close()
