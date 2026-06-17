"""
Unit tests for the frozen number-formatting contract used in the HTML reports.

``format_number_to_string`` is implemented (identically, by design) in three
places:
  * bin/create_consolidated_report.py  (module-level function)
  * bin/create_single_sample_report.py (SingleSampleHTMLReport static method)
  * bin/create_multi_sample_report.py  (MultipleSampleSummaries static method)

The documented contract is: a string containing a "." is parsed as a float and
rendered to exactly 2 d.p.; otherwise the integer string is returned unchanged.
We test all three implementations against the same table so they cannot drift
apart silently.
"""

import importlib

import pytest

# Import the three implementations.
import create_consolidated_report as ccr
import create_single_sample_report as cssr
import create_multi_sample_report as cmsr

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
    ("multi_sample_staticmethod", cmsr.MultipleSampleSummaries.format_number_to_string),
]


@pytest.mark.unit
@pytest.mark.parametrize("impl_name,fn", IMPLEMENTATIONS, ids=[n for n, _ in IMPLEMENTATIONS])
@pytest.mark.parametrize("value,expected", CASES, ids=[c[0] for c in CASES])
def test_format_number_to_string(impl_name, fn, value, expected):
    assert fn(value) == expected


@pytest.mark.unit
def test_all_three_implementations_agree():
    """Guard against the three frozen copies drifting apart."""
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
