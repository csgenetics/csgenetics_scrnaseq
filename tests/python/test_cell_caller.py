"""
Unit tests for the deterministic numeric core of bin/cell_caller.py.

The CellCaller class performs all its work in __init__ (arg parsing, h5ad
reading, plotting). We do NOT run __init__. Instead we build a bare instance
with ``object.__new__`` and set only the attributes the method under test
needs, so we can exercise the pure numeric logic in isolation:

  * ``get_cutoff``         - turns a probability-density dataframe into a log10
                             count threshold (minima / inflection / default).
  * ``get_prob_dens_data`` - KDE -> evenly sampled pdf dataframe.
  * ``assign_barcode_type``- mixed-species barcode classification rules.

The plotting methods (generate_pdf_plot / generate_barnyard_plot) are not unit
tested; they emit HTML and carry no metric arithmetic.
"""

import numpy as np
import pandas as pd
import pytest

import cell_caller as cc


def _bare_caller(minimum_count_threshold=100):
    """A CellCaller instance with __init__ bypassed and just the attrs we need."""
    obj = object.__new__(cc.CellCaller)
    obj.minimum_count_threshold = minimum_count_threshold
    return obj


# ---------------------------------------------------------------------------
# get_cutoff
# ---------------------------------------------------------------------------

@pytest.mark.unit
def test_get_cutoff_picks_minimum_between_two_peaks():
    """
    A clean bimodal density with a trough above log10(min_thresh+1) should
    return the trough location (the smallest qualifying local minimum).
    """
    caller = _bare_caller(minimum_count_threshold=100)
    # min threshold gate is log10(101) ~= 2.004. Build x from 1.5 to 5.0.
    x = np.linspace(1.5, 5.0, 200)
    # Two Gaussian humps centred at 2.4 and 4.0; trough sits ~3.2 (> 2.004).
    y = np.exp(-((x - 2.4) ** 2) / (2 * 0.15 ** 2)) + np.exp(-((x - 4.0) ** 2) / (2 * 0.20 ** 2))
    pdf_df = pd.DataFrame({"data_space": x, "evaluated": y})

    cutoff = caller.get_cutoff(pdf_df)

    # The cutoff must be the trough between the two peaks.
    assert 3.0 < cutoff < 3.5
    # And it must be one of the actual sampled data_space points.
    assert cutoff in set(pdf_df["data_space"].values)


@pytest.mark.unit
def test_get_cutoff_defaults_to_minimum_threshold_when_monotonic():
    """
    A strictly monotonic (decreasing) density has no qualifying minimum and no
    inflection point with near-zero gradient above the gate, so the cutoff must
    fall back to log10(minimum_count_threshold + 1).
    """
    caller = _bare_caller(minimum_count_threshold=100)
    x = np.linspace(0.0, 5.0, 200)
    y = np.exp(-x)  # smooth monotonic decay
    pdf_df = pd.DataFrame({"data_space": x, "evaluated": y})

    cutoff = caller.get_cutoff(pdf_df)

    assert cutoff == pytest.approx(np.log10(100 + 1))


@pytest.mark.unit
def test_get_cutoff_ignores_minima_below_gate():
    """
    A trough that sits BELOW the minimum-count gate must be ignored. Here the
    only trough is at ~1.0 (log10 scale) which is below log10(101)~=2.004, so
    with no qualifying inflection point the result defaults to the gate.
    """
    caller = _bare_caller(minimum_count_threshold=100)
    x = np.linspace(0.0, 1.8, 200)
    # Two peaks at 0.4 and 1.6 -> trough ~1.0, all below the 2.004 gate.
    y = np.exp(-((x - 0.4) ** 2) / (2 * 0.12 ** 2)) + np.exp(-((x - 1.6) ** 2) / (2 * 0.12 ** 2))
    pdf_df = pd.DataFrame({"data_space": x, "evaluated": y})

    cutoff = caller.get_cutoff(pdf_df)

    assert cutoff == pytest.approx(np.log10(100 + 1))


@pytest.mark.unit
def test_get_cutoff_picks_smallest_qualifying_minimum():
    """With multiple qualifying troughs, the smallest (leftmost) is chosen."""
    caller = _bare_caller(minimum_count_threshold=100)
    x = np.linspace(2.1, 6.0, 300)
    # Three peaks -> two troughs, both above the 2.004 gate.
    y = (
        np.exp(-((x - 2.5) ** 2) / (2 * 0.12 ** 2))
        + np.exp(-((x - 4.0) ** 2) / (2 * 0.12 ** 2))
        + np.exp(-((x - 5.5) ** 2) / (2 * 0.12 ** 2))
    )
    pdf_df = pd.DataFrame({"data_space": x, "evaluated": y})

    cutoff = caller.get_cutoff(pdf_df)

    # Smallest trough is the one between peaks 1 and 2 (~3.25), not ~4.75.
    assert 3.0 < cutoff < 3.6


# ---------------------------------------------------------------------------
# get_prob_dens_data
# ---------------------------------------------------------------------------

@pytest.mark.unit
def test_get_prob_dens_data_shape_and_range():
    caller = _bare_caller()
    rng = np.random.default_rng(0)
    counts = np.concatenate([rng.normal(2.0, 0.1, 500), rng.normal(4.0, 0.1, 500)])

    pdf_df = caller.get_prob_dens_data(counts)

    # 200 evenly spaced sample points across the data range.
    assert len(pdf_df) == 200
    assert pdf_df["data_space"].min() == pytest.approx(counts.min())
    assert pdf_df["data_space"].max() == pytest.approx(counts.max())
    # Density values are non-negative.
    assert (pdf_df["evaluated"] >= 0).all()


# ---------------------------------------------------------------------------
# assign_barcode_type (mixed species)
# ---------------------------------------------------------------------------

def _barcode_classifier(hsap_thres, mmus_thres):
    obj = object.__new__(cc.CellCaller)
    obj.hsap_thres = hsap_thres
    obj.mmus_thres = mmus_thres
    return obj


@pytest.mark.unit
@pytest.mark.parametrize(
    "hsap,mmus,expected",
    [
        (200, 200, "multiplet"),     # both strictly above
        (10, 10, "noise"),           # both strictly below
        (200, 10, "single-cell"),    # hsap above, mmus below
        (10, 200, "single-cell"),    # mmus above, hsap below
        (100, 10, "single-cell"),    # hsap exactly at thresh (not > ), mmus below -> not multiplet/noise
        (100, 100, "single-cell"),   # both exactly at thresh -> neither > nor < both -> single-cell
    ],
)
def test_assign_barcode_type(hsap, mmus, expected):
    classifier = _barcode_classifier(hsap_thres=100, mmus_thres=100)
    assert classifier.assign_barcode_type(hsap, mmus) == expected
