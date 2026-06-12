"""
Unit tests for the pure arithmetic helpers in bin/summary_statistics.py.

The SummaryStatistics class does most of its work by reading many input files
in __init__, so we do not instantiate it. Instead we test the static/pure
pieces that carry the actual metric arithmetic:

  * ``as_perc``             - fraction -> percentage
  * ``get_non_zero_sum``    - sum of non-zero entries of a 1-D array
  * the sequencing-saturation formula used in ``get_duplication_stats``

The sequencing-saturation arithmetic is inlined inside a file-reading method,
so we reproduce the exact expression here and pin it; if the formula in the
script changes, this test documents the previously-frozen behaviour. The
``as_perc`` helper it composes with is the real one from the module.
"""

import numpy as np
import pytest

import summary_statistics as ss

SS = ss.SummaryStatistics


# ---------------------------------------------------------------------------
# as_perc
# ---------------------------------------------------------------------------

@pytest.mark.unit
@pytest.mark.parametrize(
    "value,expected",
    [
        (0.0, 0.0),
        (1.0, 100.0),
        (0.5, 50.0),
        (0.1234, 12.34),
        (0.999, 99.9),
    ],
)
def test_as_perc(value, expected):
    assert SS.as_perc(value) == pytest.approx(expected)


# ---------------------------------------------------------------------------
# get_non_zero_sum
# ---------------------------------------------------------------------------

@pytest.mark.unit
def test_get_non_zero_sum_ignores_zeros():
    arr = np.array([0, 5, 0, 10, 0, 2])
    # Sum of non-zero entries == sum of all (zeros add nothing) == 17.
    assert SS.get_non_zero_sum(arr) == 17


@pytest.mark.unit
def test_get_non_zero_sum_all_zero():
    arr = np.zeros(10)
    assert SS.get_non_zero_sum(arr) == 0


@pytest.mark.unit
def test_get_non_zero_sum_with_negatives():
    # np.nonzero selects all non-zero (including negatives); they sum normally.
    arr = np.array([-3, 0, 4, 0, -1])
    assert SS.get_non_zero_sum(arr) == 0  # -3 + 4 - 1 == 0


# ---------------------------------------------------------------------------
# sequencing saturation formula
# ---------------------------------------------------------------------------

def _sequencing_saturation(reads_before, reads_after):
    """Reproduces the formula in get_duplication_stats: as_perc(1 - after/before)."""
    if reads_before != 0:
        return SS.as_perc(1 - (reads_after / reads_before))
    return 0.0


@pytest.mark.unit
def test_sequencing_saturation_typical():
    # 100 reads collapse to 25 unique -> 75% saturation.
    assert _sequencing_saturation(100, 25) == pytest.approx(75.0)


@pytest.mark.unit
def test_sequencing_saturation_no_duplication():
    # No collapse -> 0% saturation.
    assert _sequencing_saturation(100, 100) == pytest.approx(0.0)


@pytest.mark.unit
def test_sequencing_saturation_zero_reads_before():
    # Guard branch: 0 input reads -> 0.0 (avoids divide-by-zero).
    assert _sequencing_saturation(0, 0) == 0.0
