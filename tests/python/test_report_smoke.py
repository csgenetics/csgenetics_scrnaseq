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

CI runs the non-browser contract in the exact configured production container,
then persists that generated report for Chromium in a separate Python 3.11 lane
matching the declared Conda report environment.
"""

import os
from pathlib import Path
import shutil
import subprocess
import sys

import plotly.graph_objects as go
import plotly.io as pio
import pytest

try:
    from playwright.sync_api import sync_playwright
except ImportError:  # production report container intentionally has no browser tooling
    sync_playwright = None

# Runs the report generator as a subprocess and renders it in a browser.
pytestmark = pytest.mark.integration

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(REPO_ROOT, "bin"))
import create_consolidated_report  # noqa: E402

GENERATOR = os.path.join(REPO_ROOT, "bin", "create_consolidated_report.py")
TEMPLATE = os.path.join(REPO_ROOT, "templates", "consolidated_report_template.html.jinja2")
VENDOR_DIR = os.path.join(REPO_ROOT, "assets", "vendor")
FIXTURE_DIR = os.path.join(os.path.dirname(__file__), "fixtures", "report")

EXPECTED_SAMPLES = ["SAMPLE1", "SAMPLE2"]
# fixture has, per sample: 1 cell-caller + 1 qc-cascade plot, plus 1 multi-sample
# qc-cascade across all samples => 2*2 + 1 = 5 Plotly figures.
EXPECTED_PLOTS = 5
ADVERSARIAL_FRAGMENTS = [
    pytest.param(
        '<script>fetch("https://network.example.test/data")</script>',
        id="fetch",
    ),
    pytest.param(
        '<script>new Image().src="https://images.example.test/pixel"</script>',
        id="image-src",
    ),
    pytest.param(
        '<img srcset="data:image/gif;base64,AAAA 1x, https://images.example.test/x 2x">',
        id="mixed-srcset",
    ),
    pytest.param(
        '<iframe srcdoc="<script>fetch(\'https://frame.example.test\')</script>"></iframe>',
        id="iframe-srcdoc",
    ),
    pytest.param(
        '<svg><use href="https://svg.example.test/icons.svg#x"></use></svg>',
        id="svg-use",
    ),
    pytest.param(
        '<base href="https://base.example.test/"><img src="relative.png">',
        id="base-relative",
    ),
    pytest.param(
        '<form action="https://forms.example.test/submit"></form>',
        id="form-action",
    ),
    pytest.param(
        '<script nonce="csgenetics-trusted-report-script">fetch("https://network.example.test")</script>',
        id="forged-csp-nonce",
    ),
    pytest.param('<script type="text/javascript">', id="unclosed-script"),
    pytest.param("<html></html>", id="trailing-empty-html-wrapper"),
    pytest.param("<body></body><head></head>", id="misordered-wrapper"),
]


def _copy_report_fixture(work):
    for name in os.listdir(FIXTURE_DIR):
        shutil.copy(os.path.join(FIXTURE_DIR, name), work / name)


def _run_report_generator(work, *, mixed=False, capture_output=False):
    return subprocess.run(
        [
            sys.executable,
            GENERATOR,
            TEMPLATE,
            str(mixed).upper(),
            VENDOR_DIR,
            str(work / "multisample_qc_cascade.html"),
        ],
        cwd=work,
        capture_output=capture_output,
        text=True,
    )


def _assert_no_report_outputs(work):
    assert not (work / "consolidated_report.html").exists()
    assert not (work / "multisample_out.csv").exists()


@pytest.fixture(scope="module")
def report_html(tmp_path_factory):
    """Generate the consolidated report from the fixture; return the .html path."""
    work = tmp_path_factory.mktemp("report")
    _copy_report_fixture(work)

    # Same invocation Nextflow uses: <template> <mixed_species> <vendor_dir> <multi_qc_cascade>
    result = _run_report_generator(work)
    assert result.returncode == 0
    html = work / "consolidated_report.html"
    assert html.exists(), "generator did not produce consolidated_report.html"
    return html


@pytest.fixture(scope="module")
def browser():
    if sync_playwright is None:
        pytest.skip("Playwright is not installed")
    with sync_playwright() as p:
        b = p.chromium.launch()
        yield b
        b.close()


@pytest.fixture
def page(browser, report_html):
    """Load the report offline and collect JS errors and network attempts."""
    context = browser.new_context(offline=True)
    pg = context.new_page()
    errors = []
    remote_requests = []
    pg.on("console", lambda m: errors.append(m.text) if m.type == "error" else None)
    pg.on("pageerror", lambda e: errors.append(str(e)))
    pg.on(
        "request",
        lambda request: remote_requests.append(request.url)
        if request.url.startswith(("http://", "https://"))
        else None,
    )
    pg.goto(report_html.as_uri())
    pg.wait_for_timeout(1200)  # let inline Plotly.newPlot() calls settle
    pg._console_errors = errors
    pg._remote_requests = remote_requests
    yield pg
    context.close()


def test_no_js_console_errors(page):
    assert page._console_errors == [], f"JS errors in rendered report: {page._console_errors}"


def test_report_renders_without_attempting_network_access(page):
    assert page._remote_requests == [], \
        f"offline report attempted external requests: {page._remote_requests}"


def test_legacy_fragment_resources_are_removed(report_html):
    """The fixture deliberately contains the Google Fonts tag emitted by older
    Cell Caller versions. It must not survive consolidation as an active tag."""
    fixture = os.path.join(FIXTURE_DIR, "SAMPLE2_counts_pdf_with_threshold.html")
    with open(fixture, encoding="utf-8") as fh:
        assert "fonts.googleapis.com" in fh.read()

    with open(report_html, encoding="utf-8") as fh:
        html = fh.read()
    assert '<link href="https://fonts.googleapis.com' not in html
    assert '<script src="https://cdn.plot.ly' not in html


def test_report_has_restrictive_csp_and_only_trusted_scripts(page, report_html):
    with open(report_html, encoding="utf-8") as fh:
        html = fh.read()
    for directive in (
        "connect-src 'none'",
        "frame-src 'none'",
        "object-src 'none'",
        "media-src 'none'",
        "base-uri 'none'",
        "form-action 'none'",
    ):
        assert directive in html, f"report CSP is missing {directive!r}"

    script_nonces = page.eval_on_selector_all(
        "script", "els => els.map(element => element.getAttribute('nonce'))"
    )
    assert script_nonces
    assert set(script_nonces) == {"csgenetics-trusted-report-script"}


def test_dense_scattergl_fragment_renders_offline_in_consolidated_report(
    tmp_path, browser
):
    """Dense barnyards use scattergl, which must survive the strict CSP."""
    for name in os.listdir(FIXTURE_DIR):
        shutil.copy(os.path.join(FIXTURE_DIR, name), tmp_path / name)

    figure = go.Figure(
        data=go.Scattergl(x=list(range(1_001)), y=list(range(1_001)))
    )
    fragment = pio.to_html(
        figure,
        full_html=False,
        include_plotlyjs=False,
        config={"responsive": True, "displaylogo": False},
    )
    (tmp_path / "SAMPLE1_counts_pdf_with_threshold.html").write_text(
        fragment, encoding="utf-8"
    )
    subprocess.run(
        [
            sys.executable,
            GENERATOR,
            TEMPLATE,
            "FALSE",
            VENDOR_DIR,
            str(tmp_path / "multisample_qc_cascade.html"),
        ],
        cwd=tmp_path,
        check=True,
    )

    context = browser.new_context(offline=True)
    pg = context.new_page()
    errors = []
    remote_requests = []
    pg.on("console", lambda message: errors.append(message.text) if message.type == "error" else None)
    pg.on("pageerror", lambda error: errors.append(str(error)))
    pg.on(
        "request",
        lambda request: remote_requests.append(request.url)
        if request.url.startswith(("http://", "https://"))
        else None,
    )
    pg.goto((tmp_path / "consolidated_report.html").as_uri(), wait_until="load")
    pg.wait_for_selector(".gl-container canvas", state="attached", timeout=10_000)

    assert remote_requests == []
    assert errors == []
    context.close()


def test_production_container_report_renders_offline_in_host_chromium(browser):
    """CircleCI passes the report made by the configured production image here."""
    report_path_value = os.environ.get("CSGENETICS_PRODUCTION_CONTAINER_REPORT")
    if not report_path_value:
        pytest.skip("no production-container report artifact was supplied")
    report_path = Path(report_path_value)
    assert report_path.is_file(), f"missing production-container report: {report_path}"

    context = browser.new_context(offline=True)
    pg = context.new_page()
    errors = []
    remote_requests = []
    pg.on("console", lambda message: errors.append(message.text) if message.type == "error" else None)
    pg.on("pageerror", lambda error: errors.append(str(error)))
    pg.on(
        "request",
        lambda request: remote_requests.append(request.url)
        if request.url.startswith(("http://", "https://"))
        else None,
    )
    pg.goto(report_path.as_uri(), wait_until="load")
    pg.wait_for_selector(".js-plotly-plot .main-svg", state="attached", timeout=10_000)

    n_divs = pg.locator(".js-plotly-plot").count()
    n_with_svg = pg.eval_on_selector_all(
        ".js-plotly-plot",
        "els => els.filter(element => element.querySelector('.main-svg')).length",
    )
    assert n_divs >= EXPECTED_PLOTS
    assert n_with_svg == n_divs, (
        f"{n_divs - n_with_svg} of {n_divs} production-container plots "
        "drew no SVG"
    )
    assert remote_requests == []
    assert errors == []
    context.close()


@pytest.mark.parametrize(
    "case,role,filename,error_text",
    [
        (
            "missing-single",
            "QC cascade fragment for sample 'SAMPLE1'",
            "SAMPLE1.qc_cascade.html",
            "is missing",
        ),
        (
            "empty-single",
            "QC cascade fragment for sample 'SAMPLE1'",
            "SAMPLE1.qc_cascade.html",
            "is empty",
        ),
        (
            "missing-multi",
            "multi-sample QC cascade fragment",
            "multisample_qc_cascade.html",
            "is missing",
        ),
        (
            "empty-multi",
            "multi-sample QC cascade fragment",
            "multisample_qc_cascade.html",
            "is empty",
        ),
    ],
)
def test_required_qc_fragment_inputs_fail_loud_without_partial_outputs(
    tmp_path, case, role, filename, error_text
):
    _copy_report_fixture(tmp_path)
    fragment = tmp_path / filename
    if case.startswith("missing"):
        fragment.unlink()
    else:
        fragment.write_text("", encoding="utf-8")

    result = _run_report_generator(tmp_path, capture_output=True)

    assert result.returncode != 0
    assert role in result.stderr
    assert filename in result.stderr
    assert error_text in result.stderr
    _assert_no_report_outputs(tmp_path)


@pytest.mark.parametrize("optional_role", ["cell-caller", "barnyard"])
def test_existing_zero_byte_optional_plot_sentinels_are_accepted(
    tmp_path, optional_role
):
    _copy_report_fixture(tmp_path)
    if optional_role == "cell-caller":
        (tmp_path / "SAMPLE1_counts_pdf_with_threshold.html").write_text(
            "", encoding="utf-8"
        )
        mixed = False
    else:
        for sample_id in EXPECTED_SAMPLES:
            (tmp_path / f"{sample_id}_barnyard_plot.html").write_text(
                "", encoding="utf-8"
            )
        mixed = True

    result = _run_report_generator(tmp_path, mixed=mixed, capture_output=True)

    assert result.returncode == 0, result.stderr
    assert (tmp_path / "consolidated_report.html").is_file()
    assert (tmp_path / "multisample_out.csv").is_file()


@pytest.mark.parametrize(
    "mixed,filename,role",
    [
        (
            False,
            "SAMPLE1_counts_pdf_with_threshold.html",
            "Cell Caller fragment for sample 'SAMPLE1'",
        ),
        (
            True,
            "SAMPLE1_barnyard_plot.html",
            "barnyard fragment for sample 'SAMPLE1'",
        ),
    ],
)
def test_optional_plot_sentinel_must_still_exist(
    tmp_path, mixed, filename, role
):
    _copy_report_fixture(tmp_path)
    fragment = tmp_path / filename
    if fragment.exists():
        fragment.unlink()

    result = _run_report_generator(tmp_path, mixed=mixed, capture_output=True)

    assert result.returncode != 0
    assert role in result.stderr
    assert filename in result.stderr
    assert "is missing" in result.stderr
    _assert_no_report_outputs(tmp_path)


def test_non_file_fragment_fails_loud_without_partial_outputs(tmp_path):
    _copy_report_fixture(tmp_path)
    fragment = tmp_path / "SAMPLE1.qc_cascade.html"
    fragment.unlink()
    fragment.mkdir()

    result = _run_report_generator(tmp_path, capture_output=True)

    assert result.returncode != 0
    assert "QC cascade fragment for sample 'SAMPLE1'" in result.stderr
    assert fragment.name in result.stderr
    assert "not a regular file" in result.stderr
    _assert_no_report_outputs(tmp_path)


def test_unreadable_fragment_fails_loud_without_partial_outputs(tmp_path):
    _copy_report_fixture(tmp_path)
    fragment = tmp_path / "SAMPLE1.qc_cascade.html"
    fragment.chmod(0)
    try:
        result = _run_report_generator(tmp_path, capture_output=True)
    finally:
        fragment.chmod(0o600)

    assert result.returncode != 0
    assert "QC cascade fragment for sample 'SAMPLE1'" in result.stderr
    assert fragment.name in result.stderr
    assert "not readable" in result.stderr
    _assert_no_report_outputs(tmp_path)


@pytest.mark.parametrize("attack", ADVERSARIAL_FRAGMENTS)
def test_generator_rejects_unsafe_fragment_before_writing_outputs(tmp_path, attack):
    _copy_report_fixture(tmp_path)

    fragment = tmp_path / "SAMPLE1_counts_pdf_with_threshold.html"
    original = fragment.read_text(encoding="utf-8")
    fragment.write_text(
        original.replace("</body>", attack + "</body>", 1), encoding="utf-8"
    )
    result = subprocess.run(
        [
            sys.executable,
            GENERATOR,
            TEMPLATE,
            "FALSE",
            VENDOR_DIR,
            str(tmp_path / "multisample_qc_cascade.html"),
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
    )

    assert result.returncode != 0, "unsafe fragment was silently accepted"
    assert "Cell Caller fragment for sample 'SAMPLE1'" in result.stderr
    assert "is unsafe" in result.stderr
    assert fragment.name in result.stderr
    _assert_no_report_outputs(tmp_path)


def test_csp_blocks_active_content_if_a_finished_report_is_tampered(report_html, browser, tmp_path):
    tampered = tmp_path / "tampered_report.html"
    with open(report_html, encoding="utf-8") as fh:
        html = fh.read()
    attack = (
        '<script>window.__csp_script_ran=true; '
        'fetch("https://network.example.test/data")</script>'
        '<img src="https://images.example.test/pixel" '
        'onerror="window.__csp_handler_ran=true">'
    )
    tampered.write_text(html.replace("</body>", attack + "</body>"), encoding="utf-8")

    context = browser.new_context()
    pg = context.new_page()
    remote_requests = []
    blocked_requests = []
    pg.on(
        "request",
        lambda request: remote_requests.append(request.url)
        if request.url.startswith(("http://", "https://"))
        else None,
    )
    pg.on(
        "requestfailed",
        lambda request: blocked_requests.append((request.url, request.failure)),
    )
    pg.goto(tampered.as_uri())
    pg.wait_for_timeout(500)

    assert pg.evaluate("window.__csp_script_ran === undefined")
    assert pg.evaluate("window.__csp_handler_ran === undefined")
    # Chromium exposes the blocked image through its request event, but the CSP
    # cancels it before a network response. The injected script never runs, so
    # its fetch is not attempted at all.
    assert remote_requests == ["https://images.example.test/pixel"]
    assert blocked_requests == [
        ("https://images.example.test/pixel", "csp")
    ]
    context.close()


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


def test_tower_report_mappings_match_published_outputs():
    tower_path = os.path.join(REPO_ROOT, "tower.yml")
    with open(tower_path, encoding="utf-8") as fh:
        tower = fh.read()

    expected = {
        "report/consolidated_report.html": "CS Genetics scRNA-seq report",
        "report/multisample_out.csv": "Cross-sample metrics (CSV)",
        "pipeline_info/execution_report.html": "Nextflow execution report",
        "pipeline_info/execution_timeline.html": "Nextflow execution timeline",
        "pipeline_info/execution_trace.txt": "Nextflow execution trace",
        "pipeline_info/pipeline_dag.html": "Nextflow workflow diagram",
    }
    for path, display in expected.items():
        mapping = f'  "{path}":\n    display: "{display}"'
        assert mapping in tower, f"tower.yml is missing the exact mapping {path!r}"

    for stale_path in (
        "multisample_report.html",
        "*_report.html",
        "execution_timeline_*.html",
        "execution_report_*.html",
    ):
        assert stale_path not in tower, f"tower.yml still contains stale path {stale_path!r}"


def test_ci_exercises_configured_production_report_image_and_browser_handoff():
    images_config = Path(REPO_ROOT, "conf", "images.config").read_text(encoding="utf-8")
    circle_config = Path(REPO_ROOT, ".circleci", "config.yml").read_text(
        encoding="utf-8"
    )
    configured_image = "quay.io/csgenetics/html_build:0.1.0"
    assert f"container = '{configured_image}'" in images_config

    production_start = circle_config.index("  report-production-container:")
    production_end = circle_config.index(
        "  # The declared Conda report environment", production_start
    )
    production_job = circle_config[production_start:production_end]
    assert f"- image: {configured_image}" in production_job
    assert "python -m pip install --quiet pytest" in production_job
    assert 'plotly==' not in production_job
    assert "persist_to_workspace" in production_job

    assert "- image: cimg/python:3.11" in circle_config
    assert "attach_workspace" in circle_config
    assert "CSGENETICS_PRODUCTION_CONTAINER_REPORT=" in circle_config
    assert "report-smoke:\n          requires:\n            - report-production-container" in circle_config


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


def test_run_provenance_base64_treats_shell_characters_as_data(tmp_path):
    """The production Nextflow boundary must not interpolate provenance in shell."""
    import base64
    import json

    marker = tmp_path / "shell-payload-executed"
    for name in os.listdir(FIXTURE_DIR):
        shutil.copy(os.path.join(FIXTURE_DIR, name), tmp_path / name)
    run_name = f"O'Brien;$(touch {marker})"
    provenance = json.dumps({
        "genome": "GRCh38", "annotation": "customer's annotation.gtf",
        "mixed": False, "pipeline_ver": "2.0.0", "commit": "abc1234",
        "revision": "devel", "run_name": run_name, "session_id": "session-1",
        "start": "2026-08-12T12:00:00Z", "nf_version": "26.04.1",
        "outdir": "/results/O'Brien; echo inert", "barcode_kit": "kit.csv",
        "count_threshold": 100,
        "homepage": "https://github.com/csgenetics/csgenetics_scrnaseq",
    })
    encoded = base64.b64encode(provenance.encode("utf-8")).decode("ascii")

    subprocess.run(
        [sys.executable, GENERATOR, TEMPLATE, "FALSE", VENDOR_DIR,
         str(tmp_path / "multisample_qc_cascade.html"), f"base64:{encoded}"],
        cwd=tmp_path, check=True,
    )

    html = (tmp_path / "consolidated_report.html").read_text(encoding="utf-8")
    assert run_name.replace("'", "&#39;") in html
    assert "/results/O&#39;Brien; echo inert" in html
    assert not marker.exists()


def test_nextflow_report_provenance_is_encoded_before_shell_interpolation():
    module = (Path(REPO_ROOT) / "modules/local/consolidated_report/main.nf").read_text(
        encoding="utf-8"
    )
    workflow = (Path(REPO_ROOT) / "main.nf").read_text(encoding="utf-8")
    assert "val(provenance_base64)" in module
    assert "'base64:${provenance_base64}'" in module
    assert "provenance_json.getBytes('UTF-8').encodeBase64().toString()" in workflow
    assert "'${provenance_json}'" not in module


@pytest.mark.parametrize(
    "encoded", ["base64:not!base64", "base64:/w=="]
)
def test_invalid_provenance_base64_fails_loud(encoded):
    with pytest.raises(ValueError, match="run provenance base64 is invalid"):
        create_consolidated_report.load_run_provenance(
            ["generator", "template", "FALSE", "vendor", "multi", encoded]
        )


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
