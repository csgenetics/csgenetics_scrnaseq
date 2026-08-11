"""
Tests for the deterministic numeric core and CLI behavior of bin/cell_caller.py.

The CellCaller class performs all its work in __init__ (arg parsing, h5ad
reading, plotting). We do NOT run __init__. Instead we build a bare instance
with ``object.__new__`` and set only the attributes the method under test
needs, so we can exercise the pure numeric logic in isolation:

  * ``get_cutoff``         - turns a probability-density dataframe into a log10
                             count threshold (minima / inflection / default).
  * ``get_prob_dens_data`` - KDE -> evenly sampled pdf dataframe.
  * ``assign_barcode_type``- mixed-species barcode classification rules.

The integration cases execute the CLI against tiny synthetic h5ad files so the
manual-threshold value emitted to the downstream filtering process is covered.
"""

import os
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest
from scipy.sparse import csr_matrix

import cell_caller as cc


SCRIPT = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "bin", "cell_caller.py")
)


def _bare_caller(minimum_count_threshold=100):
    """A CellCaller instance with __init__ bypassed and just the attrs we need."""
    obj = object.__new__(cc.CellCaller)
    obj.minimum_count_threshold = minimum_count_threshold
    return obj


def _write_single_species_h5ad(path, total_counts):
    """Write a tiny one-gene matrix whose row sums equal total_counts."""
    obs_names = [f"barcode_{index}" for index in range(len(total_counts))]
    adata = cc.ad.AnnData(
        X=csr_matrix(np.asarray(total_counts, dtype=np.float64).reshape(-1, 1)),
        obs=pd.DataFrame(index=obs_names),
        var=pd.DataFrame(index=["gene_1"]),
    )
    adata.write_h5ad(path)


def _write_mixed_species_h5ad(path):
    """Write populations that are valid for filtering but singular for KDE."""
    obs = pd.DataFrame(
        {
            "hsap_counts": [5, 5, 1, 1],
            "mmus_counts": [1, 1, 7, 7],
        },
        index=["human_1", "human_2", "mouse_1", "mouse_2"],
    )
    adata = cc.ad.AnnData(
        X=csr_matrix(np.ones((4, 1), dtype=np.float64)),
        obs=obs,
        var=pd.DataFrame(index=["gene_1"]),
    )
    adata.write_h5ad(path)


def _run_cell_caller(tmp_path, count_matrix, manual_threshold, single_species=True):
    return subprocess.run(
        [
            sys.executable,
            SCRIPT,
            "--sample_name",
            "sample",
            "--single_species",
            str(single_species),
            "--minimum_count_threshold",
            "100",
            "--count_matrix",
            str(count_matrix),
            "--manual_threshold_str",
            manual_threshold,
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
    )


# ---------------------------------------------------------------------------
# manual-threshold precedence and validation
# ---------------------------------------------------------------------------


@pytest.mark.integration
@pytest.mark.parametrize("total_counts", [[7, 7, 7], [7]])
def test_manual_threshold_is_used_for_degenerate_and_small_inputs(tmp_path, total_counts):
    count_matrix = tmp_path / "counts.h5ad"
    _write_single_species_h5ad(count_matrix, total_counts)
    manual_log_threshold = str(np.log10(11))

    result = _run_cell_caller(tmp_path, count_matrix, manual_log_threshold)

    assert result.returncode == 0, result.stderr
    assert result.stdout == "10"
    assert (tmp_path / "sample_counts_pdf_with_threshold.html").stat().st_size == 0
    assert (tmp_path / "sample_pdf_with_cutoff.html").stat().st_size == 0


@pytest.mark.integration
def test_zero_manual_threshold_is_valid_and_authoritative(tmp_path):
    count_matrix = tmp_path / "counts.h5ad"
    _write_single_species_h5ad(count_matrix, [7])

    result = _run_cell_caller(tmp_path, count_matrix, "0")

    assert result.returncode == 0, result.stderr
    assert result.stdout == "0"


@pytest.mark.integration
def test_manual_mixed_thresholds_survive_degenerate_species_distributions(tmp_path):
    count_matrix = tmp_path / "counts.h5ad"
    _write_mixed_species_h5ad(count_matrix)
    hsap_log_threshold = np.log10(11)
    mmus_log_threshold = np.log10(21)

    result = _run_cell_caller(
        tmp_path,
        count_matrix,
        f"{hsap_log_threshold}_{mmus_log_threshold}",
        single_species=False,
    )

    assert result.returncode == 0, result.stderr
    assert result.stdout == "10_20"
    assert (tmp_path / "sample_hsap_pdf_with_cutoff.html").stat().st_size == 0
    assert (tmp_path / "sample_mmus_pdf_with_cutoff.html").stat().st_size == 0
    assert (tmp_path / "sample_barnyard_plot.html").stat().st_size > 0


@pytest.mark.integration
def test_partial_manual_mixed_threshold_preserves_manual_and_automatic_fallback(tmp_path):
    count_matrix = tmp_path / "counts.h5ad"
    _write_mixed_species_h5ad(count_matrix)

    result = _run_cell_caller(
        tmp_path,
        count_matrix,
        f"{np.log10(11)}_nan",
        single_species=False,
    )

    assert result.returncode == 0, result.stderr
    assert result.stdout == "10_100"


@pytest.mark.integration
def test_degenerate_automatic_threshold_keeps_minimum_fallback(tmp_path):
    count_matrix = tmp_path / "counts.h5ad"
    _write_single_species_h5ad(count_matrix, [7])

    result = _run_cell_caller(tmp_path, count_matrix, "nan")

    assert result.returncode == 0, result.stderr
    assert result.stdout == "100"


@pytest.mark.integration
@pytest.mark.parametrize("manual_threshold", ["-0.1", "inf", "NaN", "invalid", "309"])
def test_invalid_manual_threshold_fails_clearly(tmp_path, manual_threshold):
    result = _run_cell_caller(
        tmp_path,
        tmp_path / "not-read-because-validation-fails.h5ad",
        manual_threshold,
    )

    assert result.returncode == 2
    assert "manual Cell Caller threshold" in result.stderr
    assert result.stdout == ""


@pytest.mark.unit
def test_malformed_mixed_manual_threshold_fails_clearly(monkeypatch, capsys):
    caller = object.__new__(cc.CellCaller)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            SCRIPT,
            "--sample_name",
            "sample",
            "--single_species",
            "false",
            "--count_matrix",
            "unused.h5ad",
            "--manual_threshold_str",
            "2.5",
        ],
    )

    with pytest.raises(SystemExit) as error:
        caller.parse_arguments()

    assert error.value.code == 2
    assert "exactly two values separated by an underscore" in capsys.readouterr().err


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
