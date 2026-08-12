#!/usr/bin/env python

"""
Create ONE consolidated, self-contained HTML report for a CS Genetics
scRNA-seq experiment.

This replaces the previous per-sample reports (``${sample_id}_report.html``)
and the multi-sample report family (``multisample_report.html``,
``multisample_summary_plots.html``) with a single navigable, styled report:

  * an experiment overview header (sample count + headline metrics),
  * a per-sample view (selectable via Bootstrap tabs) reproducing everything
    the old per-sample report showed: headline cards, the Read QC /
    Deduplication / Cell metrics groups, the cell-caller plot, the qc-cascade
    plot, and (mixed-species) the barnyard plot,
  * the experiment-wide multi-sample summary plots (violin/box) and the
    experiment-wide qc-cascade plot,
  * the cross-sample metrics table (same content as multisample_out.csv).

HARD CONSTRAINT: this script does NOT change any metric VALUES. It reads the
SAME ``${sample_id}.metrics.csv`` files produced by summary_statistics.py and
renders them. Number formatting matches the old reports exactly
(floats -> 2 d.p., ints unchanged). It also re-emits ``multisample_out.csv``
with the identical schema and values as before.

OFFLINE-SAFE: Plotly.js, Bootstrap CSS/JS, Bootstrap-icons and the Lexend font
are embedded inline from locally vendored assets. The rendered report contains
no external CDN URLs and renders with no network.

Usage:
  create_consolidated_report.py <template> <mixed_species> <vendor_dir> \\
      <multi_qc_cascade_html> [provenance_json]

Nextflow supplies provenance in the optional positional argument with a
``base64:`` prefix. Groovy performs the encoding before task generation, so raw
customer paths and run names are never interpolated into shell source. Direct
JSON remains supported for backwards-compatible invocation.

All per-sample inputs are staged FLAT into the working directory by Nextflow.
Samples are discovered from the ``*.metrics.csv`` files, and each sample's plot
fragments are located by their conventional filenames in the same directory:
  - <sample_id>.metrics.csv
  - <sample_id>_counts_pdf_with_threshold.html   (cell-caller fragment)
  - <sample_id>_barnyard_plot.html               (barnyard fragment; mixed only)
  - <sample_id>.qc_cascade.html                  (qc-cascade fragment)
"""

import sys
import base64
import binascii
import os
import json
import re
import stat
from html.parser import HTMLParser
from pathlib import Path
from collections import defaultdict
from urllib.parse import urlsplit

from jinja2 import Environment
from plotly.offline import get_plotlyjs

from create_single_sample_report import get_cell_stat_cat_dict_obj


# Seqera Platform's documented report limits are expressed in decimal MB.  Keep
# these values here (rather than in the Nextflow wrapper) so the warning is based
# on the file that was actually written.
SEQERA_PREVIEW_LIMIT_BYTES = 10_000_000
SEQERA_DOWNLOAD_LIMIT_BYTES = 25_000_000


# Plot fragments are inserted with Jinja's ``safe`` filter, so accepting arbitrary
# HTML here would turn an offline report into an execution/network boundary.  The
# contract below accepts only the exact structural subset emitted by
# ``plotly.io.to_html(..., full_html=False, include_plotlyjs=False)``. Plotly
# upgrades must revalidate these expressions and the pinned offline browser
# suite before changing the production dependency. The only
# legacy markup removed is a document wrapper, the old Google Fonts stylesheet,
# and Plotly's own CDN loader. Everything else fails loud.
REPORT_SCRIPT_NONCE = "csgenetics-trusted-report-script"
_PLOT_DIV_ID_RE = re.compile(r"[A-Za-z0-9_-]+\Z")
_PLOTLY_SCRIPT_PREFIX_RE = re.compile(
    r"\A\s*window\.PLOTLYENV\s*=\s*window\.PLOTLYENV\s*\|\|\s*\{\}\s*;\s*"
    r"if\s*\(\s*document\.getElementById\s*\(\s*"
)
_PLOTLY_SCRIPT_MIDDLE_RE = re.compile(
    r"\s*\)\s*\)\s*\{\s*Plotly\.newPlot\s*\(\s*"
)
_PLOTLY_SCRIPT_SUFFIX_RE = re.compile(r"\s*\)\s*;?\s*\}\s*;?\s*\Z")
_REMOTE_OR_EXECUTABLE_URL_RE = re.compile(
    r"(?:https?:)?//|\b(?:javascript|vbscript|file|ftp|wss?):",
    re.IGNORECASE,
)
_RESOURCE_JSON_KEYS = {
    "href",
    "images",
    "plotlyserverurl",
    "source",
    "src",
    "srcset",
    "topojsonurl",
    "url",
}
# Every trace constructor used by this pipeline. Plotly Express automatically
# selects WebGL-backed ``scattergl`` for barnyard plots above 1,000 barcodes.
_ALLOWED_TRACE_TYPES = {"bar", "box", "scatter", "scattergl"}


def format_number_to_string(number_str):
    """
    Format numbers exactly as the previous reports did:
    floats -> 2 d.p.; ints unchanged. Frozen formatting contract.
    Used ONLY for the multisample_out.csv data file (must stay byte-identical;
    thousands separators would corrupt the comma-delimited file).
    """
    if "." in number_str:
        return f"{float(number_str):.2f}"
    else:
        return number_str


def format_number_for_display(number_str):
    """
    Human-readable DISPLAY formatting for the on-screen tables and headline cards:
    add thousands separators to the integer part and keep the existing decimals
    (floats stay at 2 d.p., matching the value contract -- no rounding, so
    percentages/rates/fractions are unchanged; e.g. 181004808 -> "181,004,808",
    12106.54 -> "12,106.54", 94.30 -> "94.30"). DISPLAY ONLY -- never used for the
    CSV. Non-numeric sentinels (e.g. "nan", "nan_nan") pass through unchanged.
    """
    s = number_str.strip()
    try:
        if "." in s:
            return f"{float(s):,.2f}"
        return f"{int(s):,}"
    except ValueError:
        return number_str


class PlotFragmentError(ValueError):
    """A nonempty report fragment violated the trusted Plotly contract."""


def _normalise_network_url(url):
    return f"https:{url}" if url.startswith("//") else url


def _is_legacy_google_font(attrs):
    """Return True only for the Google Fonts link emitted by older releases."""
    if set(attrs) - {"href", "rel", "type"}:
        return False
    href = attrs.get("href")
    rel = attrs.get("rel", "")
    if not href or "stylesheet" not in rel.lower().split():
        return False
    parsed = urlsplit(_normalise_network_url(href))
    return (
        parsed.scheme in {"http", "https"}
        and parsed.hostname == "fonts.googleapis.com"
        and parsed.path in {"/css", "/css2"}
    )


def _is_legacy_plotly_cdn_script(attrs):
    """Return True only for Plotly's generated external library loader."""
    if set(attrs) - {"charset", "crossorigin", "integrity", "src", "type"}:
        return False
    src = attrs.get("src")
    if not src:
        return False
    parsed = urlsplit(_normalise_network_url(src))
    return (
        parsed.scheme in {"http", "https"}
        and parsed.hostname == "cdn.plot.ly"
        and re.fullmatch(r"/plotly(?:-[A-Za-z0-9.]+)?\.min\.js", parsed.path) is not None
    )


def _validate_json_safety(value, context):
    """Reject Plotly JSON fields capable of loading or executing external data."""
    if isinstance(value, dict):
        for key, child in value.items():
            key_lower = str(key).lower()
            if key_lower in _RESOURCE_JSON_KEYS or key_lower.endswith("src"):
                raise PlotFragmentError(
                    f"{context} contains forbidden resource-bearing key {key!r}"
                )
            _validate_json_safety(child, f"{context}.{key}")
    elif isinstance(value, list):
        for index, child in enumerate(value):
            _validate_json_safety(child, f"{context}[{index}]")
    elif isinstance(value, str) and _REMOTE_OR_EXECUTABLE_URL_RE.search(value):
        raise PlotFragmentError(f"{context} contains a remote or executable URL")


def _decode_json_at(script, position, label):
    position += len(script[position:]) - len(script[position:].lstrip())
    try:
        return json.JSONDecoder().raw_decode(script, position)
    except json.JSONDecodeError as exc:
        raise PlotFragmentError(
            f"Plotly.newPlot {label} is not strict JSON ({exc.msg})"
        ) from exc


def _validate_plotly_script(script):
    """Validate one exact Plotly ``to_html`` initializer and return its div id."""
    prefix = _PLOTLY_SCRIPT_PREFIX_RE.match(script)
    if prefix is None:
        raise PlotFragmentError(
            "inline script is not a generated Plotly.newPlot initializer"
        )

    lookup_id, position = _decode_json_at(script, prefix.end(), "lookup id")
    if not isinstance(lookup_id, str):
        raise PlotFragmentError("Plotly lookup id must be a JSON string")

    middle = _PLOTLY_SCRIPT_MIDDLE_RE.match(script, position)
    if middle is None:
        raise PlotFragmentError("Plotly initializer has an unexpected preamble")
    position = middle.end()

    labels = ("target id", "data", "layout", "config")
    values = []
    for index, label in enumerate(labels):
        value, position = _decode_json_at(script, position, label)
        values.append(value)
        if index < len(labels) - 1:
            comma = re.match(r"\s*,\s*", script[position:])
            if comma is None:
                raise PlotFragmentError(
                    f"Plotly.newPlot {label} is not followed by the expected comma"
                )
            position += comma.end()

    if _PLOTLY_SCRIPT_SUFFIX_RE.match(script, position) is None:
        raise PlotFragmentError("Plotly initializer contains trailing executable code")

    target_id, data, layout, config = values
    if not isinstance(target_id, str) or target_id != lookup_id:
        raise PlotFragmentError("Plotly lookup id and target id do not match")
    if not isinstance(data, list) or not data:
        raise PlotFragmentError("Plotly data must be a nonempty JSON array")
    if not isinstance(layout, dict) or not isinstance(config, dict):
        raise PlotFragmentError("Plotly layout and config must be JSON objects")

    for index, trace in enumerate(data):
        if not isinstance(trace, dict):
            raise PlotFragmentError(f"Plotly trace {index} is not a JSON object")
        trace_type = trace.get("type")
        if trace_type not in _ALLOWED_TRACE_TYPES:
            raise PlotFragmentError(
                f"Plotly trace {index} has unsupported offline type {trace_type!r}"
            )

    _validate_json_safety(data, "Plotly data")
    _validate_json_safety(layout, "Plotly layout")
    _validate_json_safety(config, "Plotly config")
    return target_id


class _TrustedPlotlyFragmentParser(HTMLParser):
    """Strict structural parser that preserves accepted source bytes.

    ``HTMLParser`` identifies token boundaries; edits are applied to those exact
    boundaries after validation. Script bodies are never searched/replaced as
    HTML, so strings containing markup cannot be silently rewritten.
    """

    def __init__(self, fragment, script_nonce):
        super().__init__(convert_charrefs=False)
        self.fragment = fragment
        self.script_nonce = script_nonce
        self.line_offsets = [0]
        self.line_offsets.extend(match.end() for match in re.finditer("\n", fragment))
        self.stack = []
        self.edits = []
        self.plot_ids = set()
        self.script_targets = set()
        self.wrapper_counts = {"html": 0, "head": 0, "body": 0}
        self.doctype_seen = False
        self.wrapper_mode = None
        self.wrapper_phase = "start"

    def _offset(self):
        line, column = self.getpos()
        return self.line_offsets[line - 1] + column

    def _fail(self, message):
        line, column = self.getpos()
        raise PlotFragmentError(f"line {line}, column {column}: {message}")

    def _attribute_map(self, attrs):
        names = [name for name, _value in attrs]
        if len(names) != len(set(names)):
            self._fail("duplicate HTML attribute")
        if any(value is None for _name, value in attrs):
            self._fail("boolean HTML attributes are not valid in a Plotly fragment")
        return dict(attrs)

    def _record_wrapper(self, tag, attrs, start, end):
        if tag == "html":
            if set(attrs) - {"lang"}:
                self._fail("legacy <html> wrapper has unexpected attributes")
            if self.stack or self.wrapper_mode is not None:
                self._fail("<html> wrapper must be the outermost element")
            self.wrapper_mode = "html"
            self.wrapper_phase = "before-head-or-body"
        elif attrs:
            self._fail(f"legacy <{tag}> wrapper has unexpected attributes")

        if tag == "head":
            if not self.stack:
                if self.doctype_seen or self.wrapper_mode is not None:
                    self._fail(
                        "top-level <head> must begin the legacy head-then-body wrapper"
                    )
                self.wrapper_mode = "head-body"
                self.wrapper_phase = "head-open"
            elif (
                len(self.stack) == 1
                and self.stack[0]["tag"] == "html"
                and self.wrapper_mode == "html"
                and self.wrapper_phase == "before-head-or-body"
            ):
                self.wrapper_phase = "head-open"
            else:
                self._fail("<head> must precede <body> directly inside <html>")
        elif tag == "body":
            if (
                len(self.stack) == 1
                and self.stack[0]["tag"] == "html"
                and self.wrapper_mode == "html"
                and self.wrapper_phase in {"before-head-or-body", "after-head"}
            ):
                self.wrapper_phase = "body-open"
            elif (
                not self.stack
                and self.wrapper_mode == "head-body"
                and self.wrapper_phase == "after-head"
            ):
                self.wrapper_phase = "body-open"
            else:
                self._fail(
                    "<body> must follow optional <head> inside <html>, or a top-level <head>"
                )

        self.wrapper_counts[tag] += 1
        if self.wrapper_counts[tag] > 1:
            self._fail(f"multiple <{tag}> wrappers are not valid")
        self.edits.append((start, end, ""))
        self.stack.append({"tag": tag, "kind": "wrapper", "start": start, "start_end": end})

    def _record_content_start(self, tag):
        """Require plot content to be bare or inside the one supported body."""
        if not self.stack:
            if self.doctype_seen and self.wrapper_mode is None:
                self._fail("<!doctype html> must be followed by an outer <html> wrapper")
            if self.wrapper_mode is None:
                self.wrapper_mode = "bare"
                self.wrapper_phase = "content"
            elif self.wrapper_mode != "bare":
                self._fail(f"element <{tag}> appears outside the completed document wrapper")
            return

        wrapper_parent = next(
            (item["tag"] for item in reversed(self.stack) if item["kind"] == "wrapper"),
            None,
        )
        if wrapper_parent == "html":
            self._fail(f"element <{tag}> must be inside the document <body>")
        if wrapper_parent == "head" and tag not in {"link", "script"}:
            self._fail(f"element <{tag}> is not a permitted legacy <head> resource")

    def _validate_plot_style(self, style):
        """Accept only Plotly's inert inline height/width declarations."""
        if style is None:
            return

        style_parts = {}
        for declaration in style.split(";"):
            if not declaration.strip():
                continue
            if ":" not in declaration:
                self._fail("Plotly div style is malformed")
            key, value = (part.strip().lower() for part in declaration.split(":", 1))
            if key in style_parts:
                self._fail(f"duplicate Plotly div style property {key!r}")
            style_parts[key] = value

        if not style_parts or set(style_parts) - {"height", "width"} or any(
            re.fullmatch(r"(?:100%|[0-9]+(?:\.[0-9]+)?px)", value) is None
            for value in style_parts.values()
        ):
            self._fail("Plotly div style contains an unsafe declaration")

    def _record_div(self, attrs, start, end):
        if self.stack and self.stack[-1]["kind"] == "plot":
            self._fail("Plotly graph div must be empty")
        if self.stack and self.stack[-1]["tag"] == "head":
            self._fail("plot div cannot appear inside <head>")

        # Plotly 5.x emitted a bare outer container while newer releases may put
        # the same inert height/width style used by the graph div on that outer
        # container.  Both are generated forms; class/id remain exclusive to the
        # actual graph landmark below.
        if set(attrs) <= {"style"}:
            self._validate_plot_style(attrs.get("style"))
            kind = "container"
        else:
            if set(attrs) - {"class", "id", "style"}:
                self._fail("Plotly graph div has unexpected attributes")
            classes = attrs.get("class", "").split()
            plot_id = attrs.get("id", "")
            if classes != ["plotly-graph-div"] or not _PLOT_DIV_ID_RE.fullmatch(plot_id):
                self._fail("non-container div is not a Plotly graph landmark")
            if plot_id in self.plot_ids:
                self._fail(f"duplicate Plotly graph id {plot_id!r}")

            self._validate_plot_style(attrs.get("style"))

            self.plot_ids.add(plot_id)
            kind = "plot"

        self.stack.append({"tag": "div", "kind": kind, "start": start, "start_end": end})

    def handle_starttag(self, tag, attrs):
        start = self._offset()
        raw = self.get_starttag_text()
        end = start + len(raw)
        attrs = self._attribute_map(attrs)

        if self.stack and self.stack[-1]["kind"] == "plot":
            self._fail("Plotly graph div must be empty")

        if tag in self.wrapper_counts:
            self._record_wrapper(tag, attrs, start, end)
            return

        self._record_content_start(tag)

        if tag == "link":
            if not (self.stack and self.stack[-1]["tag"] == "head"):
                self._fail("legacy stylesheet link must be inside <head>")
            if not _is_legacy_google_font(attrs):
                self._fail("external or unknown <link> is forbidden")
            self.edits.append((start, end, ""))
            return

        if tag == "div":
            self._record_div(attrs, start, end)
            return

        if tag == "script":
            src = attrs.get("src")
            if src is not None:
                if not _is_legacy_plotly_cdn_script(attrs):
                    self._fail("external or unknown script loader is forbidden")
                kind = "legacy-script"
            else:
                if set(attrs) - {"type"} or attrs.get("type", "text/javascript") not in {
                    "application/javascript",
                    "text/javascript",
                }:
                    self._fail("inline Plotly script has unexpected attributes")
                if self.stack and self.stack[-1]["tag"] == "head":
                    self._fail("inline fragment script cannot appear inside <head>")
                kind = "plot-script"
            self.stack.append({"tag": tag, "kind": kind, "start": start, "start_end": end})
            return

        self._fail(f"element <{tag}> is forbidden by the Plotly fragment contract")

    def handle_startendtag(self, tag, attrs):
        if tag != "link":
            self._fail(f"self-closing <{tag}> is forbidden by the Plotly fragment contract")
        self.handle_starttag(tag, attrs)

    def handle_endtag(self, tag):
        start = self._offset()
        end = self.fragment.find(">", start)
        if end == -1:
            self._fail(f"unterminated </{tag}> tag")
        end += 1
        raw = self.fragment[start:end]
        if re.fullmatch(rf"</{re.escape(tag)}\s*>", raw, re.IGNORECASE) is None:
            self._fail(f"malformed </{tag}> tag")
        if not self.stack or self.stack[-1]["tag"] != tag:
            expected = self.stack[-1]["tag"] if self.stack else "no open element"
            self._fail(f"mismatched </{tag}>; expected {expected}")

        opened = self.stack.pop()
        if opened["kind"] == "wrapper":
            if tag == "head":
                if self.wrapper_phase != "head-open":
                    self._fail("legacy <head> closed in an invalid wrapper state")
                self.wrapper_phase = "after-head"
            elif tag == "body":
                if self.wrapper_phase != "body-open":
                    self._fail("legacy <body> closed in an invalid wrapper state")
                self.wrapper_phase = (
                    "after-body" if self.wrapper_mode == "html" else "closed"
                )
            elif tag == "html":
                if self.wrapper_phase != "after-body":
                    self._fail("outer <html> must contain a completed <body>")
                self.wrapper_phase = "closed"
            self.edits.append((start, end, ""))
        elif opened["kind"] == "legacy-script":
            body = self.fragment[opened["start_end"]:start]
            if body.strip():
                self._fail("legacy Plotly CDN loader contains unexpected inline code")
            self.edits.append((opened["start"], end, ""))
        elif opened["kind"] == "plot-script":
            body = self.fragment[opened["start_end"]:start]
            target = _validate_plotly_script(body)
            if target in self.script_targets:
                self._fail(f"duplicate Plotly initializer for graph {target!r}")
            self.script_targets.add(target)
            self.edits.append(
                (opened["start_end"] - 1, opened["start_end"] - 1,
                 f' nonce="{self.script_nonce}"')
            )

    def handle_data(self, data):
        if self.stack and self.stack[-1]["kind"] in {"plot-script", "legacy-script"}:
            return
        if data.strip():
            self._fail("non-whitespace text is forbidden outside Plotly scripts")

    def handle_entityref(self, name):
        self._fail(f"entity reference &{name}; is forbidden outside Plotly JSON")

    def handle_charref(self, name):
        self._fail(f"character reference &#{name}; is forbidden outside Plotly JSON")

    def handle_comment(self, data):
        self._fail("HTML comments are forbidden in Plotly fragments")

    def handle_decl(self, decl):
        if (
            self.doctype_seen
            or decl.strip().lower() != "doctype html"
            or self.stack
            or self.wrapper_mode is not None
        ):
            self._fail("only one outer <!doctype html> declaration is allowed")
        start = self._offset()
        end = self.fragment.find(">", start)
        if end == -1:
            self._fail("unterminated doctype")
        self.edits.append((start, end + 1, ""))
        self.doctype_seen = True

    def handle_pi(self, data):
        self._fail("processing instructions are forbidden in Plotly fragments")

    def unknown_decl(self, data):
        self._fail("unknown declarations are forbidden in Plotly fragments")

    def validate_and_render(self):
        try:
            self.feed(self.fragment)
            self.close()
        except PlotFragmentError:
            raise
        except Exception as exc:
            raise PlotFragmentError(f"HTML parser rejected fragment: {exc}") from exc

        if self.stack:
            unclosed = ", ".join(f"<{item['tag']}>" for item in self.stack)
            raise PlotFragmentError(f"unclosed element(s): {unclosed}")
        if self.wrapper_mode == "html" and self.wrapper_phase != "closed":
            raise PlotFragmentError(
                "outer <html> wrapper must contain optional <head> followed by <body>"
            )
        if self.wrapper_mode == "head-body" and self.wrapper_phase != "closed":
            raise PlotFragmentError(
                "top-level legacy wrapper must contain <head> followed by <body>"
            )
        if self.doctype_seen and self.wrapper_mode != "html":
            raise PlotFragmentError(
                "<!doctype html> is only valid immediately before an outer <html> wrapper"
            )
        if not self.plot_ids:
            raise PlotFragmentError("fragment contains no Plotly graph div")
        if self.script_targets != self.plot_ids:
            missing = sorted(self.plot_ids - self.script_targets)
            unknown = sorted(self.script_targets - self.plot_ids)
            raise PlotFragmentError(
                f"Plotly graph/initializer mismatch (missing={missing}, unknown={unknown})"
            )

        previous_end = -1
        for start, end, _replacement in sorted(self.edits):
            if start < previous_end:
                raise PlotFragmentError("internal error: overlapping structural edits")
            previous_end = max(previous_end, end)

        rendered = self.fragment
        for start, end, replacement in sorted(self.edits, reverse=True):
            rendered = rendered[:start] + replacement + rendered[end:]
        return rendered


def validate_plot_fragment(fragment, script_nonce=REPORT_SCRIPT_NONCE):
    """Validate and minimally unwrap one trusted generated Plotly fragment."""
    return _TrustedPlotlyFragmentParser(fragment, script_nonce).validate_and_render()


def seqera_report_size_warning(report_size_bytes):
    """Describe how a large report can be accessed in Seqera Platform.

    Returns ``None`` while the report is previewable.  The complete report is
    always published; the separately mapped ``multisample_out.csv`` is the
    lightweight in-Platform fallback when the HTML is too large to preview.
    """
    if report_size_bytes < SEQERA_PREVIEW_LIMIT_BYTES:
        return None

    size_mb = report_size_bytes / 1_000_000
    if report_size_bytes < SEQERA_DOWNLOAD_LIMIT_BYTES:
        action = "download the HTML from the Reports tab and open it locally"
    else:
        action = "retrieve the HTML from its published output path and open it locally"

    return (
        f"WARNING: consolidated_report.html is {size_mb:.1f} MB. "
        "Seqera Platform previews only reports smaller than 10 MB and directly "
        "downloads only reports smaller than 25 MB; "
        f"{action}. Cross-sample metrics remain available in the separately "
        "mapped multisample_out.csv report."
    )


def read_fragment(
    path,
    *,
    role,
    script_nonce=REPORT_SCRIPT_NONCE,
    allow_empty=False,
):
    """Read and validate one required report input.

    Every declared input must exist, be a readable regular file, and contain a
    trusted Plotly fragment. ``allow_empty`` is reserved for the Cell Caller and
    barnyard zero-byte sentinels; it never makes a missing path optional.
    """
    display_path = "<not provided>" if path is None else os.fspath(path)
    input_label = f"{role} at {display_path!r}"
    if path is None:
        raise PlotFragmentError(f"{input_label} is missing")

    path_obj = Path(path)
    try:
        file_stat = path_obj.stat()
    except FileNotFoundError as exc:
        raise PlotFragmentError(f"{input_label} is missing") from exc
    except OSError as exc:
        raise PlotFragmentError(
            f"{input_label} cannot be inspected: {exc.strerror or exc}"
        ) from exc

    if not stat.S_ISREG(file_stat.st_mode):
        raise PlotFragmentError(f"{input_label} is not a regular file")
    read_bits = stat.S_IRUSR | stat.S_IRGRP | stat.S_IROTH
    if not file_stat.st_mode & read_bits or not os.access(path_obj, os.R_OK):
        raise PlotFragmentError(f"{input_label} is not readable")
    if file_stat.st_size == 0:
        if allow_empty:
            return None
        raise PlotFragmentError(f"{input_label} is empty; a Plotly fragment is required")

    try:
        with path_obj.open("r", encoding="utf-8") as f:
            fragment = f.read()
    except (OSError, UnicodeError) as exc:
        raise PlotFragmentError(f"{input_label} cannot be read: {exc}") from exc

    # Guard the stat/open race: a producer truncating the file between those
    # operations is still an empty required input, not a missing plot.
    if not fragment:
        if allow_empty:
            return None
        raise PlotFragmentError(f"{input_label} is empty; a Plotly fragment is required")
    try:
        return validate_plot_fragment(fragment, script_nonce)
    except PlotFragmentError as exc:
        raise PlotFragmentError(f"{input_label} is unsafe: {exc}") from exc


def read_text_asset(vendor_dir, name):
    with open(os.path.join(vendor_dir, name), "r", encoding="utf-8") as f:
        return f.read()


def load_run_provenance(argv):
    """Strictly decode and parse optional run metadata."""
    raw = argv[5] if len(argv) > 5 else ""
    if not raw or not raw.strip():
        return {}
    if raw.startswith("base64:"):
        try:
            raw = base64.b64decode(raw.removeprefix("base64:"), validate=True).decode(
                "utf-8", errors="strict"
            )
        except (binascii.Error, UnicodeDecodeError) as exc:
            raise ValueError(f"run provenance base64 is invalid: {exc}") from exc
    try:
        provenance = json.loads(raw)
    except json.JSONDecodeError as exc:
        raise ValueError(f"run provenance is not valid JSON: {exc}") from exc
    if not isinstance(provenance, dict):
        raise ValueError("run provenance must be a JSON object")
    return provenance


class ConsolidatedReport:
    def __init__(self):
        self.template_path = sys.argv[1]
        self.mixed = sys.argv[2].upper() == "TRUE"
        self.vendor_dir = sys.argv[3]
        self.multi_qc_cascade_path = sys.argv[4]
        # Nextflow passes optional run-provenance JSON as a shell-inert base64
        # argument. The direct JSON form remains supported. Absent -> no section.
        self.provenance = load_run_provenance(sys.argv)
        # All per-sample inputs are staged flat into the cwd.
        self.work_dir = "."
        self.script_nonce = REPORT_SCRIPT_NONCE

        # autoescape=True so all data fields (metric values, human-readable names,
        # tooltips, and especially the customer-supplied sample ids) are HTML-escaped.
        # The trusted inlined assets and Plotly plot fragments are raw HTML and carry
        # an explicit `| safe` in the template; nothing else is trusted.
        with open(self.template_path) as fh:
            env = Environment(autoescape=True)
            self.jinja_template = env.from_string(fh.read())

        # Discover samples by their metrics csv. Sort for a stable, logical order.
        self.sample_ids = self._discover_sample_ids()

        # Build per-sample data structures (cards + plot fragments).
        self.samples = [self._build_sample(sid) for sid in self.sample_ids]
        self.multi_qc_cascade_fragment = read_fragment(
            self.multi_qc_cascade_path,
            role="multi-sample QC cascade fragment",
            script_nonce=self.script_nonce,
        )

        # Build experiment-wide structures (cross-sample table + summary plots),
        # re-emitting multisample_out.csv with the unchanged schema/values. All
        # fragments have already passed validation before any output is written.
        self.multi_metrics_dict, self.cell_metric_tooltip_dict = self._build_multi_metrics_and_csv()
        self.summary_plot_fragment = self._build_summary_plot_fragment()

        # Headline experiment metrics for the overview header.
        self.overview = self._build_overview()

        # Vendored, inlined assets (offline-safe).
        self.assets = self._load_assets()

        self.render_and_write()

    # ------------------------------------------------------------------ #
    # Sample discovery / per-sample data
    # ------------------------------------------------------------------ #
    def _metrics_csv_for(self, sample_id):
        return os.path.join(self.work_dir, f"{sample_id}.metrics.csv")

    def _discover_sample_ids(self):
        # Discover from the flat-staged metrics csvs. Sort for a stable order.
        csvs = sorted(Path(self.work_dir).glob("*.metrics.csv"))
        return [c.name[: -len(".metrics.csv")] for c in csvs]

    def _build_sample(self, sample_id):
        # Per-sample nested metrics dict, mirroring create_single_sample_report.py:
        # metrics_dict[classification][variable_name] = (human_readable, value, tooltip)
        metrics_dict = defaultdict(dict)
        with open(self._metrics_csv_for(sample_id)) as fh:
            next(fh)  # header
            for line in fh:
                var_name, var_value, human, tooltip, group = line.strip().split(",")
                metrics_dict[group][var_name] = (human, format_number_for_display(var_value), tooltip)

        cell_plot = read_fragment(
            os.path.join(self.work_dir, f"{sample_id}_counts_pdf_with_threshold.html"),
            role=f"Cell Caller fragment for sample {sample_id!r}",
            script_nonce=self.script_nonce,
            allow_empty=True,
        )
        qc_cascade = read_fragment(
            os.path.join(self.work_dir, f"{sample_id}.qc_cascade.html"),
            role=f"QC cascade fragment for sample {sample_id!r}",
            script_nonce=self.script_nonce,
        )
        if self.mixed:
            barnyard = read_fragment(
                os.path.join(self.work_dir, f"{sample_id}_barnyard_plot.html"),
                role=f"barnyard fragment for sample {sample_id!r}",
                script_nonce=self.script_nonce,
                allow_empty=True,
            )
        else:
            barnyard = None

        return {
            "sample_id": sample_id,
            # plain dict so Jinja's `in` / attribute access behaves predictably
            "metrics_dict": {k: v for k, v in metrics_dict.items()},
            "cell_plot": cell_plot,
            "show_cell_plot": cell_plot is not None,
            "qc_cascade": qc_cascade,
            "show_qc_cascade": qc_cascade is not None,
            "barnyard": barnyard,
            "show_barnyard": barnyard is not None,
        }

    # ------------------------------------------------------------------ #
    # Experiment-wide data (reuses the frozen multi-sample logic)
    # ------------------------------------------------------------------ #
    def _build_multi_metrics_and_csv(self):
        """
        Parse the per-sample metrics csvs (no recomputation) to build both the
        re-emitted multisample_out.csv (byte-identical schema/values) and the
        cross-sample table. The CSV keeps the frozen format_number_to_string
        contract; the table gets thousands separators (see the return).
        """
        # metrics_dict[classification][(variable_name, human, description, classification)]
        #   = [value_per_sample, ...] in sample order
        metrics_dict = defaultdict(lambda: defaultdict(list))

        # Build from the SAME csvs, in the SAME sorted order as the table columns.
        csv_paths = [self._metrics_csv_for(sid) for sid in self.sample_ids]

        # Store the RAW per-sample values; format per-purpose below. The CSV keeps
        # the frozen format_number_to_string contract; the on-screen table gets
        # thousands separators. The two must not be conflated -- separators would
        # corrupt the comma-delimited CSV.
        for csv_path in csv_paths:
            with open(csv_path) as fh:
                next(fh)  # header
                for line in fh:
                    variable_name, value, human, description, classification = line.strip().split(",")
                    metrics_dict[classification][(variable_name, human, description, classification)].append(value)

        # Re-emit multisample_out.csv with the unchanged schema and values
        # (format_number_to_string: floats -> 2 d.p., ints unchanged; NO separators).
        with open("multisample_out.csv", "w") as csv_out:
            csv_out.write("variable_name,human_readable_name,description,classification,")
            csv_out.write(",".join(self.sample_ids) + "\n")
            for classification, subdict in metrics_dict.items():
                for (variable_name, human, description, _cls), vals in subdict.items():
                    csv_out.write(f"{','.join([variable_name, human, description, classification])},")
                    csv_out.write(",".join(format_number_to_string(v) for v in vals) + "\n")

        # Cell-metric tooltips for the cross-sample table (mixed only).
        cell_tooltip = get_cell_stat_cat_dict_obj(self.mixed)
        cell_tooltip = {k: v[0] for k, v in cell_tooltip.items()}

        # The cross-sample TABLE shows display-formatted values (thousands
        # separators); the CSV above kept the raw frozen contract.
        table_dict = {
            cls: {key: [format_number_for_display(v) for v in vals] for key, vals in sub.items()}
            for cls, sub in metrics_dict.items()
        }
        return table_dict, cell_tooltip

    def _build_summary_plot_fragment(self):
        """
        DISABLED (2026-06, by decision): the across-sample "Key metric distributions"
        violin/box summary plot is intentionally not generated. A violin/KDE drawn
        over the small number of samples these reports carry implies a continuous
        distribution that does not exist. Returning None hides the section cleanly
        via the template's `{% if summary_plot %}` guard (and the template block has
        also been removed). Reinstate by restoring this method from git history if a
        suitable small-N representation is ever chosen.
        """
        return None

    # ------------------------------------------------------------------ #
    # Overview header
    # ------------------------------------------------------------------ #
    def _build_overview(self):
        """A small set of headline experiment metrics for the top of the report."""
        # Per-sample headline metric keys (same as the old per-sample headline cards).
        if self.mixed:
            # In mixed-species mode summary_statistics.py classifies each cell metric under its own
            # name (e.g. num_cells), NOT under "Cell metrics" (which only holds counts_in/out_of_cells).
            headline_keys = [
                ("num_cells_total", "num_cells"),
                ("raw_reads_per_cell_total", "raw_reads_per_cell"),
                ("median_genes_detected_per_cell_total", "median_genes_detected_per_cell"),
            ]
        else:
            headline_keys = [
                ("num_cells", "Cell metrics"),
                ("raw_reads_per_cell", "Cell metrics"),
                ("median_genes_detected_per_cell", "Cell metrics"),
            ]
        return {"n_samples": len(self.sample_ids), "headline_keys": headline_keys}

    # ------------------------------------------------------------------ #
    # Vendored assets (inlined for offline-safety)
    # ------------------------------------------------------------------ #
    def _load_assets(self):
        plotlyjs = get_plotlyjs()  # full Plotly.js source, embedded once
        return {
            "bootstrap_css": read_text_asset(self.vendor_dir, "bootstrap.min.css"),
            "bootstrap_icons_css": read_text_asset(self.vendor_dir, "bootstrap-icons.embedded.css"),
            "lexend_css": read_text_asset(self.vendor_dir, "lexend.embedded.css"),
            "bootstrap_js": read_text_asset(self.vendor_dir, "bootstrap.bundle.min.js"),
            "plotly_js": plotlyjs,
            "logo_data_uri": read_text_asset(self.vendor_dir, "cs_logo_data_uri.txt").strip(),
        }

    # ------------------------------------------------------------------ #
    # Render
    # ------------------------------------------------------------------ #
    def render_and_write(self):
        cell_stat_cat_dict = get_cell_stat_cat_dict_obj(self.mixed)

        html = self.jinja_template.render(
            mixed=self.mixed,
            assets=self.assets,
            overview=self.overview,
            samples=self.samples,
            sample_name_list=self.sample_ids,
            multi_metrics_dict=self.multi_metrics_dict,
            cell_metric_tooltip_dict=self.cell_metric_tooltip_dict,
            cell_stat_cat_dict=cell_stat_cat_dict,
            summary_plot=self.summary_plot_fragment,
            multi_qc_cascade=self.multi_qc_cascade_fragment,
            provenance=self.provenance,
            csp_nonce=self.script_nonce,
        )

        with open("consolidated_report.html", "w") as fh:
            fh.write(html)

        warning = seqera_report_size_warning(os.path.getsize("consolidated_report.html"))
        if warning:
            print(warning, file=sys.stderr)


if __name__ == "__main__":
    ConsolidatedReport()
