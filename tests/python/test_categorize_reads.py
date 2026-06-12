"""
Unit tests for bin/categorize_reads.py.

These exercise the pure, importable helpers that do the in-cell vs out-of-cell
categorisation arithmetic and the barcode parsing, using tiny in-memory
anndata objects and crafted read names. The BAM-iterating part of ``main()``
is not unit-tested here (it needs a real BAM); the count-based arithmetic it
relies on (``get_counts_by_cell_status``) and the cell-extraction logic
(``get_cells_from_h5ad``) are tested directly.
"""

import json

import numpy as np
import pandas as pd
import pytest
from scipy.sparse import csr_matrix
import anndata as ad

import categorize_reads as cr


# ---------------------------------------------------------------------------
# extract_barcode_from_read_name
# ---------------------------------------------------------------------------

@pytest.mark.unit
def test_extract_barcode_basic():
    read_name = "VH01429:149:22253TGNX:1:2601:69579:24995_TGAGCCACATCGA_"
    # 13 bp barcode between the underscores
    assert cr.extract_barcode_from_read_name(read_name, 13) == "TGAGCCACATCGA"


@pytest.mark.unit
def test_extract_barcode_with_n_base():
    # The pattern allows N as a valid barcode base.
    read_name = "instr:1:flow:1:1:1:1_ACGTNACGTNACG_"
    assert cr.extract_barcode_from_read_name(read_name, 13) == "ACGTNACGTNACG"


@pytest.mark.unit
def test_extract_barcode_wrong_length_returns_none():
    # Barcode in the name is 13 bp but we ask for 12 -> the anchored pattern
    # cannot match a 12-mer flanked by underscores, so None.
    read_name = "instr:1:flow:1:1:1:1_TGAGCCACATCGA_"
    assert cr.extract_barcode_from_read_name(read_name, 12) is None


@pytest.mark.unit
def test_extract_barcode_no_barcode_returns_none():
    read_name = "instr:1:flow:1:1:1:1_no_barcode_here"
    assert cr.extract_barcode_from_read_name(read_name, 13) is None


# ---------------------------------------------------------------------------
# get_raw_reads_from_fastp
# ---------------------------------------------------------------------------

@pytest.mark.unit
def test_get_raw_reads_from_fastp(tmp_path):
    fastp = tmp_path / "fastp.json"
    fastp.write_text(json.dumps({"summary": {"before_filtering": {"total_reads": 12345}}}))
    assert cr.get_raw_reads_from_fastp(str(fastp)) == 12345


# ---------------------------------------------------------------------------
# helpers to build a tiny anndata
# ---------------------------------------------------------------------------

def _make_adata(barcodes, counts_matrix, obs_cols=None):
    """counts_matrix: list-of-lists (n_barcodes x n_genes)."""
    X = csr_matrix(np.array(counts_matrix, dtype=np.float32))
    obs = pd.DataFrame(index=barcodes)
    if obs_cols:
        for col, vals in obs_cols.items():
            obs[col] = vals
    var = pd.DataFrame(index=[f"gene{i}" for i in range(X.shape[1])])
    return ad.AnnData(X=X, obs=obs, var=var)


# ---------------------------------------------------------------------------
# get_cells_from_h5ad
# ---------------------------------------------------------------------------

@pytest.mark.unit
def test_get_cells_single_species():
    adata = _make_adata(
        ["S1_AAA", "S1_BBB", "S1_CCC"],
        [[1, 2], [3, 4], [5, 6]],
        obs_cols={"is_single_cell": [True, False, True]},
    )
    cells = cr.get_cells_from_h5ad(adata, mixed_species=False)
    assert cells["total"] == {"S1_AAA", "S1_CCC"}
    assert "hsap" not in cells


@pytest.mark.unit
def test_get_cells_mixed_species():
    adata = _make_adata(
        ["S1_AAA", "S1_BBB", "S1_CCC", "S1_DDD"],
        [[1, 2], [3, 4], [5, 6], [7, 8]],
        obs_cols={
            "is_single_cell": [True, True, False, False],
            "is_hsap_cell": [True, False, False, False],
            "is_mmus_cell": [False, True, False, False],
        },
    )
    cells = cr.get_cells_from_h5ad(adata, mixed_species=True)
    assert cells["total"] == {"S1_AAA", "S1_BBB"}
    assert cells["hsap"] == {"S1_AAA"}
    assert cells["mmus"] == {"S1_BBB"}


@pytest.mark.unit
def test_get_cells_missing_column_raises():
    adata = _make_adata(["S1_AAA"], [[1, 2]])
    with pytest.raises(ValueError, match="is_single_cell"):
        cr.get_cells_from_h5ad(adata, mixed_species=False)


# ---------------------------------------------------------------------------
# get_counts_by_cell_status
# ---------------------------------------------------------------------------

@pytest.mark.unit
def test_get_counts_by_cell_status_arithmetic():
    # Barcode totals: 3, 7, 11, 15. Cells are rows 0 and 2 -> in = 3 + 11 = 14;
    # out = 7 + 15 = 22.
    adata = _make_adata(
        ["S1_AAA", "S1_BBB", "S1_CCC", "S1_DDD"],
        [[1, 2], [3, 4], [5, 6], [7, 8]],
        obs_cols={"is_single_cell": [True, False, True, False]},
    )
    in_cells, out_cells = cr.get_counts_by_cell_status(adata, "S1")
    assert in_cells == 14
    assert out_cells == 22
    # Return type is plain int.
    assert isinstance(in_cells, int) and isinstance(out_cells, int)


@pytest.mark.unit
def test_get_counts_all_cells():
    adata = _make_adata(
        ["S1_AAA", "S1_BBB"],
        [[10, 0], [0, 5]],
        obs_cols={"is_single_cell": [True, True]},
    )
    in_cells, out_cells = cr.get_counts_by_cell_status(adata, "S1")
    assert in_cells == 15
    assert out_cells == 0


@pytest.mark.unit
def test_get_counts_missing_column_raises():
    adata = _make_adata(["S1_AAA"], [[1, 2]])
    with pytest.raises(ValueError, match="is_single_cell"):
        cr.get_counts_by_cell_status(adata, "S1")
