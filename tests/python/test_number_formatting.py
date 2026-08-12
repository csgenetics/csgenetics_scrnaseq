"""
Unit tests for the frozen number-formatting contract used in the HTML reports.

``format_number_to_string`` is implemented (identically, by design) in the two
report scripts the pipeline uses:
  * bin/create_consolidated_report.py  (module-level function)
  * bin/create_single_sample_report.py (SingleSampleHTMLReport static method)

The documented contract is: a string containing a "." is parsed as a float and
rendered to exactly 2 d.p.; otherwise the integer string is returned unchanged.
We test both implementations against the same table so they cannot drift apart
silently.
"""

import re
from pathlib import Path

import pytest

# Both report formatters are lean (jinja2 + base64 only), always importable.
import create_consolidated_report as ccr
import create_single_sample_report as cssr

# (input, expected) pairs covering the contract.
CASES = [
    ("100", "100"),          # plain int unchanged
    ("0", "0"),
    ("1234567", "1234567"),
    ("3.14159", "3.14"),     # float rounded to 2 d.p.
    ("3.1", "3.10"),         # float padded to 2 d.p.
    ("0.0", "0.00"),
    ("2.005", "2.00"),       # banker's-ish rounding via format spec
    ("99.999", "100.00"),
    ("-5", "-5"),
    ("-2.5", "-2.50"),
]

IMPLEMENTATIONS = [
    ("consolidated_module_fn", ccr.format_number_to_string),
    ("single_sample_staticmethod", cssr.SingleSampleHTMLReport.format_number_to_string),
]


@pytest.mark.unit
@pytest.mark.parametrize("impl_name,fn", IMPLEMENTATIONS, ids=[n for n, _ in IMPLEMENTATIONS])
@pytest.mark.parametrize("value,expected", CASES, ids=[c[0] for c in CASES])
def test_format_number_to_string(impl_name, fn, value, expected):
    assert fn(value) == expected


@pytest.mark.unit
def test_all_implementations_agree():
    """Guard against the two frozen formatter copies drifting apart."""
    for value, _ in CASES:
        results = {fn(value) for _, fn in IMPLEMENTATIONS}
        assert len(results) == 1, f"implementations disagree on {value!r}: {results}"


# ---------------------------------------------------------------------------
# Display formatting (HTML tables/cards only): thousands separators on the
# integer part, decimals preserved, no rounding. Kept SEPARATE from the CSV
# contract above (commas would corrupt the comma-delimited multisample_out.csv).
# ---------------------------------------------------------------------------
DISPLAY_CASES = [
    ("181004808", "181,004,808"),   # large int gets separators
    ("100", "100"),                 # small int unchanged
    ("0", "0"),
    ("-5000", "-5,000"),
    ("12106.54", "12,106.54"),      # float: separators on integer part, decimals kept
    ("94.30", "94.30"),             # percentage-like float: value unchanged, 2 d.p.
    ("3.1", "3.10"),                # padded to 2 d.p. like the frozen contract
    ("0.0", "0.00"),
    ("nan", "nan"),                 # non-numeric sentinel passes through
    ("nan_nan", "nan_nan"),         # mixed-species sentinel passes through
]


@pytest.mark.unit
@pytest.mark.parametrize("value,expected", DISPLAY_CASES, ids=[c[0] for c in DISPLAY_CASES])
def test_format_number_for_display(value, expected):
    assert ccr.format_number_for_display(value) == expected


@pytest.mark.unit
def test_display_and_csv_formatters_diverge_only_on_separators():
    """The CSV must stay separator-free (machine-readable); the display adds commas."""
    assert ccr.format_number_to_string("181004808") == "181004808"      # CSV: no separators
    assert ccr.format_number_for_display("181004808") == "181,004,808"  # display: separators


# ---------------------------------------------------------------------------
# Offline fragment handling and Seqera size guidance.
# ---------------------------------------------------------------------------


FIXTURE_FRAGMENT = (
    Path(__file__).parent / "fixtures" / "report" / "SAMPLE1_counts_pdf_with_threshold.html"
)
BARE_FIXTURE_FRAGMENT = (
    Path(__file__).parent / "fixtures" / "report" / "SAMPLE1.qc_cascade.html"
)
SCRIPT_BODY_RE = re.compile(r"<script\b[^>]*>(.*?)</script\s*>", re.IGNORECASE | re.DOTALL)


@pytest.mark.unit
def test_trusted_fragment_contract_only_removes_known_legacy_wrappers():
    raw = FIXTURE_FRAGMENT.read_text(encoding="utf-8")
    # Plotly's historical CDN loader is the second and only other legacy network
    # element the contract permits. Its body must remain empty.
    raw = raw.replace(
        "</head>",
        '<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script></head>',
        1,
    )

    cleaned = ccr.validate_plot_fragment(raw)

    assert "fonts.googleapis.com" not in cleaned
    assert "cdn.plot.ly" not in cleaned
    assert "<head" not in cleaned
    assert "<body" not in cleaned
    assert f'nonce="{ccr.REPORT_SCRIPT_NONCE}"' in cleaned

    raw_inline_scripts = [
        body for body in SCRIPT_BODY_RE.findall(raw) if "Plotly.newPlot" in body
    ]
    cleaned_inline_scripts = SCRIPT_BODY_RE.findall(cleaned)
    assert cleaned_inline_scripts == raw_inline_scripts, \
        "structural validation must not rewrite valid Plotly script bytes"


@pytest.mark.unit
@pytest.mark.parametrize(
    "wrapped",
    [
        lambda bare: bare,
        lambda bare: f"<head></head><body>{bare}</body>",
        lambda bare: f"<html><body>{bare}</body></html>",
        lambda bare: (
            "<!doctype html><html lang=\"en\"><head>"
            '<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Lexend">'
            f"</head><body>{bare}</body></html>"
        ),
    ],
    ids=["bare", "legacy-head-body", "outer-html", "doctype-html-head-body"],
)
def test_trusted_fragment_contract_accepts_only_supported_wrapper_shapes(wrapped):
    bare = BARE_FIXTURE_FRAGMENT.read_text(encoding="utf-8")
    trusted = ccr.validate_plot_fragment(wrapped(bare))
    assert "Plotly.newPlot" in trusted
    assert "<html" not in trusted.lower()
    assert "<head" not in trusted.lower()
    assert "<body" not in trusted.lower()


@pytest.mark.unit
@pytest.mark.parametrize(
    "wrapped",
    [
        lambda bare: bare + "<html></html>",
        lambda bare: bare + "<head></head>",
        lambda bare: bare + "<body></body>",
        lambda bare: bare + "<head></head><body></body>",
        lambda bare: f"<body>{bare}</body><head></head>",
        lambda bare: f"<head></head><body><head></head>{bare}</body>",
        lambda bare: f"<html><body>{bare}</body><head></head></html>",
        lambda bare: f"<html><head></head>{bare}<body></body></html>",
        lambda bare: "<!doctype html>" + bare,
        lambda bare: f"<html>{bare}</html>",
        lambda bare: f"<head></head>{bare}<body></body>",
        lambda bare: f"<html><body>{bare}</body></html><body></body>",
    ],
    ids=[
        "trailing-empty-html",
        "trailing-empty-head",
        "trailing-empty-body",
        "trailing-head-body",
        "body-before-head",
        "head-nested-in-body",
        "head-after-body-in-html",
        "plot-outside-body-in-html",
        "doctype-without-html",
        "html-without-body",
        "plot-between-top-level-head-body",
        "wrapper-after-closed-html",
    ],
)
def test_trusted_fragment_contract_rejects_misordered_or_trailing_wrappers(wrapped):
    bare = BARE_FIXTURE_FRAGMENT.read_text(encoding="utf-8")
    with pytest.raises(ccr.PlotFragmentError):
        ccr.validate_plot_fragment(wrapped(bare))


@pytest.mark.unit
def test_plotly_json_markup_like_string_is_preserved_byte_for_byte():
    bare = BARE_FIXTURE_FRAGMENT.read_text(encoding="utf-8")
    marker = "literal <body><!-- still JSON text --></body>"
    raw = bare.replace('"name":"Input reads"', f'"name":"{marker}"', 1)

    trusted = ccr.validate_plot_fragment(raw)

    assert marker in trusted
    raw_body = SCRIPT_BODY_RE.findall(raw)
    trusted_body = SCRIPT_BODY_RE.findall(trusted)
    assert trusted_body == raw_body


@pytest.mark.unit
@pytest.mark.parametrize(
    "attack",
    [
        '<script>fetch("https://network.example.test/data")</script>',
        '<script>new Image().src="https://images.example.test/pixel"</script>',
        '<img srcset="data:image/gif;base64,AAAA 1x, https://images.example.test/x 2x">',
        '<iframe srcdoc="<script>fetch(\'https://frame.example.test\')</script>"></iframe>',
        '<svg><use href="https://svg.example.test/icons.svg#x"></use></svg>',
        '<base href="https://base.example.test/"><img src="relative.png">',
        '<form action="https://forms.example.test/submit"></form>',
        f'<script nonce="{ccr.REPORT_SCRIPT_NONCE}">fetch("https://network.example.test")</script>',
        '<script type="text/javascript">',
    ],
    ids=[
        "fetch",
        "image-src",
        "mixed-srcset",
        "iframe-srcdoc",
        "svg-use",
        "base-relative",
        "form-action",
        "forged-nonce",
        "unclosed-script",
    ],
)
def test_trusted_fragment_contract_rejects_active_or_malformed_html(attack):
    raw = FIXTURE_FRAGMENT.read_text(encoding="utf-8")
    attacked = raw.replace("</body>", attack + "</body>", 1)
    with pytest.raises(ccr.PlotFragmentError):
        ccr.validate_plot_fragment(attacked)


@pytest.mark.unit
def test_trusted_fragment_contract_rejects_remote_plotly_json():
    raw = FIXTURE_FRAGMENT.read_text(encoding="utf-8")
    attacked = raw.replace(
        '"name":""',
        '"name":"https://network.example.test/from-plot-data"',
        1,
    )
    with pytest.raises(ccr.PlotFragmentError, match="remote or executable URL"):
        ccr.validate_plot_fragment(attacked)


@pytest.mark.unit
def test_read_fragment_failure_names_the_bad_path(tmp_path):
    fragment_path = tmp_path / "broken_plot.html"
    fragment_path.write_text("<script>", encoding="utf-8")

    with pytest.raises(ccr.PlotFragmentError, match="broken_plot.html"):
        ccr.read_fragment(fragment_path, role="test plot fragment")


@pytest.mark.unit
@pytest.mark.parametrize(
    "size,expected_fragment",
    [
        (ccr.SEQERA_PREVIEW_LIMIT_BYTES - 1, None),
        (ccr.SEQERA_PREVIEW_LIMIT_BYTES, "download the HTML from the Reports tab"),
        (ccr.SEQERA_DOWNLOAD_LIMIT_BYTES - 1, "download the HTML from the Reports tab"),
        (ccr.SEQERA_DOWNLOAD_LIMIT_BYTES, "retrieve the HTML from its published output path"),
    ],
)
def test_seqera_report_size_warning(size, expected_fragment):
    warning = ccr.seqera_report_size_warning(size)
    if expected_fragment is None:
        assert warning is None
    else:
        assert expected_fragment in warning
        assert "multisample_out.csv" in warning
