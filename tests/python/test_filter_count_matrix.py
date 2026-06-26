"""
Integration test for bin/filter_count_matrix.py.

filter_count_matrix.py reads a raw count-matrix h5ad plus a count threshold and:
  * single species: flags barcodes with total_counts >= threshold as
    is_single_cell, and writes a filtered h5ad containing only those.
  * mixed species (threshold "<hsap>_<mmus>", argv[4] == "TRUE"): applies the
    hsap/mmus classification rules to set is_called_cell / is_single_cell /
    is_hsap_cell / is_mmus_cell.

The script reads sys.argv directly at import, so we drive it as a subprocess on
a crafted input h5ad and assert the obs annotations on the output.
"""

import os
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest

anndata = pytest.importorskip("anndata")
from scipy.sparse import csr_matrix  # noqa: E402

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
SCRIPT = os.path.join(REPO_ROOT, "bin", "filter_count_matrix.py")


def _make_var(n_genes):
    return pd.DataFrame(
        {
            "gene_id": [f"ENSG_{i}" for i in range(n_genes)],
            "gene_name": [f"Gene{i}" for i in range(n_genes)],
        },
        index=[f"Gene{i}" for i in range(n_genes)],
    )


@pytest.mark.integration
def test_filter_single_species(tmp_path):
    # 4 barcodes with total_counts 50, 100, 150, 99. Threshold 100 -> cells are
    # those with total_counts >= 100 -> barcodes 1 and 2.
    barcodes = ["S1_A", "S1_B", "S1_C", "S1_D"]
    X = csr_matrix(np.array(
        [[25, 25], [60, 40], [100, 50], [50, 49]], dtype=np.float32
    ))
    obs = pd.DataFrame({"total_counts": [50, 100, 150, 99]}, index=barcodes)
    adata = anndata.AnnData(X=X, obs=obs, var=_make_var(2))

    in_h5ad = tmp_path / "raw.h5ad"
    adata.write(str(in_h5ad))

    result = subprocess.run(
        [sys.executable, SCRIPT, "100", str(in_h5ad), "S1", "FALSE"],
        capture_output=True, text=True, cwd=str(tmp_path),
    )
    assert result.returncode == 0, f"STDERR:{result.stderr}"

    raw_out = tmp_path / "S1.100.raw_feature_bc_matrix.h5ad"
    filt_out = tmp_path / "S1.100.filtered_feature_bc_matrix.h5ad"
    assert raw_out.exists() and filt_out.exists()

    raw = anndata.read_h5ad(str(raw_out))
    is_cell = dict(zip(raw.obs_names, raw.obs["is_single_cell"]))
    assert is_cell["S1_A"] == False  # 50 < 100  # noqa: E712
    assert is_cell["S1_B"] == True   # 100 >= 100
    assert is_cell["S1_C"] == True   # 150 >= 100
    assert is_cell["S1_D"] == False  # 99 < 100

    filt = anndata.read_h5ad(str(filt_out))
    assert set(filt.obs_names) == {"S1_B", "S1_C"}

    # Tripartite outputs are written too.
    for fn in ["matrix.mtx.gz", "barcodes.tsv.gz", "features.tsv.gz"]:
        assert (tmp_path / fn).exists()


@pytest.mark.integration
def test_filter_mixed_species(tmp_path):
    """
    Mixed-species classification with hsap_thresh=100, mmus_thresh=100.

    Barcodes:
      H : hsap=200, mmus=10   -> hsap single-cell
      M : hsap=10,  mmus=200  -> mmus single-cell
      MULT : hsap=200, mmus=200 -> called (either >=) but NOT single-cell
      NOISE: hsap=10, mmus=10  -> not called, not single-cell
    """
    barcodes = ["S1_H", "S1_M", "S1_MULT", "S1_NOISE"]
    # Build a 2-gene matrix; the per-gene values don't matter for obs-based
    # filtering, only the hsap_counts / mmus_counts obs columns do.
    X = csr_matrix(np.ones((4, 2), dtype=np.float32))
    obs = pd.DataFrame(
        {
            "hsap_counts": [200, 10, 200, 10],
            "mmus_counts": [10, 200, 200, 10],
            "total_counts": [210, 210, 400, 20],
        },
        index=barcodes,
    )
    adata = anndata.AnnData(X=X, obs=obs, var=_make_var(2))

    in_h5ad = tmp_path / "raw_mixed.h5ad"
    adata.write(str(in_h5ad))

    result = subprocess.run(
        [sys.executable, SCRIPT, "100_100", str(in_h5ad), "S1", "TRUE"],
        capture_output=True, text=True, cwd=str(tmp_path),
    )
    assert result.returncode == 0, f"STDERR:{result.stderr}"

    raw = anndata.read_h5ad(str(tmp_path / "S1.100_100.raw_feature_bc_matrix.h5ad"))

    called = dict(zip(raw.obs_names, raw.obs["is_called_cell"]))
    single = dict(zip(raw.obs_names, raw.obs["is_single_cell"]))
    hsap = dict(zip(raw.obs_names, raw.obs["is_hsap_cell"]))
    mmus = dict(zip(raw.obs_names, raw.obs["is_mmus_cell"]))

    # H: hsap cell.
    assert hsap["S1_H"] == True and mmus["S1_H"] == False  # noqa: E712
    assert single["S1_H"] == True and called["S1_H"] == True

    # M: mmus cell.
    assert hsap["S1_M"] == False and mmus["S1_M"] == True
    assert single["S1_M"] == True and called["S1_M"] == True

    # MULT: called (>= a threshold) but both above -> not single, not hsap/mmus.
    assert called["S1_MULT"] == True
    assert single["S1_MULT"] == False
    assert hsap["S1_MULT"] == False and mmus["S1_MULT"] == False

    # NOISE: nothing.
    assert called["S1_NOISE"] == False
    assert single["S1_NOISE"] == False

    # Filtered h5ad keeps only single cells (H and M).
    filt = anndata.read_h5ad(str(tmp_path / "S1.100_100.filtered_feature_bc_matrix.h5ad"))
    assert set(filt.obs_names) == {"S1_H", "S1_M"}


@pytest.mark.integration
def test_filter_empty_input(tmp_path):
    """An empty (zero-byte) input h5ad -> empty outputs and clean exit."""
    empty = tmp_path / "empty.h5ad"
    empty.write_text("")

    result = subprocess.run(
        [sys.executable, SCRIPT, "100", str(empty), "S1", "FALSE"],
        capture_output=True, text=True, cwd=str(tmp_path),
    )
    assert result.returncode == 0, f"STDERR:{result.stderr}"
    assert (tmp_path / "S1.100.filtered_feature_bc_matrix.empty.h5ad").exists()
    assert (tmp_path / "S1.100.raw_feature_bc_matrix.empty.h5ad").exists()
