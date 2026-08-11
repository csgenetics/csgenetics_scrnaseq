"""
Unit tests for bin/categorize_reads.py.

These exercise the pure, importable helpers that do the in-cell vs out-of-cell
categorisation arithmetic and the barcode parsing, using tiny in-memory
anndata objects and crafted read names. Small synthetic BAMs also exercise
primary-alignment handling and the pipeline's intentional header-only empty
BAM contract.
"""

import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest
from scipy.sparse import csr_matrix
import anndata as ad

import categorize_reads as cr


SCRIPT = Path(__file__).resolve().parents[2] / "bin" / "categorize_reads.py"


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


def _write_bam(path, records=(), *, include_sequence_dictionary=True):
    header = {"HD": {"VN": "1.6"}}
    if include_sequence_dictionary:
        header["SQ"] = [{"SN": "chr1", "LN": 1000}]
    with cr.pysam.AlignmentFile(path, "wb", header=header) as bam:
        for query_name, flag, start in records:
            read = cr.pysam.AlignedSegment()
            read.query_name = query_name
            read.query_sequence = "ACGT"
            read.flag = flag
            read.reference_id = 0
            read.reference_start = start
            read.mapping_quality = 60
            read.cigar = ((0, 4),)
            read.query_qualities = cr.pysam.qualitystring_to_array("IIII")
            bam.write(read)
    return path


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


@pytest.mark.unit
@pytest.mark.parametrize("invalid_value", [np.nan, np.inf, -1, 0.5])
def test_get_counts_rejects_invalid_raw_count_values(invalid_value):
    adata = _make_adata(
        ["S1_AAA"],
        [[invalid_value]],
        obs_cols={"is_single_cell": [True]},
    )
    with pytest.raises(ValueError, match="count matrix"):
        cr.get_counts_by_cell_status(adata, "S1")


@pytest.mark.unit
@pytest.mark.parametrize("raw_reads", [-1, 1.5, np.nan, np.inf, True, "3"])
def test_raw_read_count_must_be_finite_nonnegative_integer(tmp_path, raw_reads):
    fastp = tmp_path / "fastp.json"
    fastp.write_text(
        json.dumps({"summary": {"before_filtering": {"total_reads": raw_reads}}})
    )
    with pytest.raises(ValueError, match="total_reads"):
        cr.get_raw_reads_from_fastp(fastp)


@pytest.mark.unit
def test_only_named_zero_byte_h5ad_is_accepted_as_empty(tmp_path):
    sentinel = tmp_path / "S1.raw_feature_bc_matrix.empty.h5ad"
    sentinel.touch()
    assert cr.load_count_matrix(sentinel) is None

    generic = tmp_path / "S1.raw_feature_bc_matrix.h5ad"
    generic.touch()
    with pytest.raises(ValueError, match=r"\*\.empty\.h5ad"):
        cr.load_count_matrix(generic)


@pytest.mark.unit
def test_named_empty_h5ad_must_be_zero_bytes(tmp_path):
    sentinel = tmp_path / "S1.raw_feature_bc_matrix.empty.h5ad"
    sentinel.write_bytes(b"not actually an empty sentinel")
    with pytest.raises(ValueError, match="must be zero bytes"):
        cr.load_count_matrix(sentinel)


@pytest.mark.unit
def test_corrupt_nonempty_h5ad_is_not_treated_as_empty(tmp_path):
    corrupt = tmp_path / "S1.raw_feature_bc_matrix.h5ad"
    corrupt.write_bytes(b"not an HDF5 file")
    with pytest.raises((OSError, ValueError)):
        cr.load_count_matrix(corrupt)


@pytest.mark.integration
def test_primary_multimapper_is_counted_once(tmp_path):
    cell_barcode = "AAAAAAAAAAAAA"
    noncell_barcode = "CCCCCCCCCCCCC"
    bam = _write_bam(
        tmp_path / "aligned.bam",
        [
            (f"read1_{cell_barcode}_", 0, 10),
            (f"read1_{cell_barcode}_", 256, 20),
            (f"read1_{cell_barcode}_", 2048, 30),
            (f"read2_{noncell_barcode}_", 0, 40),
        ],
    )

    metrics = cr.categorize_primary_reads(bam, {cell_barcode}, 13)

    assert metrics == {
        "reads_in_cells": 1,
        "reads_out_of_cells": 1,
        "total_reads_processed": 2,
        "skipped_secondary": 1,
        "skipped_supplementary": 1,
    }


@pytest.mark.integration
def test_named_empty_h5ad_and_header_only_bam_emit_finite_zero_metrics(tmp_path):
    bam = _write_bam(
        tmp_path / "S1_Aligned.sortedByCoord.out.bam",
        include_sequence_dictionary=False,
    )
    sentinel = tmp_path / "S1.100.raw_feature_bc_matrix.empty.h5ad"
    sentinel.touch()
    fastp = tmp_path / "S1.R1.preQC.fastp.json"
    fastp.write_text(
        json.dumps({"summary": {"before_filtering": {"total_reads": 0}}}),
        encoding="utf-8",
    )

    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--sample_id",
            "S1",
            "--star_bam",
            str(bam),
            "--raw_count_matrix_h5ad",
            str(sentinel),
            "--fastp_json",
            str(fastp),
            "--barcode_length",
            "13",
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stderr
    metrics = pd.read_csv(tmp_path / "S1.read_categorization.csv").iloc[0]
    for name in (
        "raw_reads",
        "reads_in_cells",
        "reads_out_of_cells",
        "unusable_reads",
        "counts_in_cells",
        "counts_out_of_cells",
    ):
        assert np.isfinite(metrics[name])
        assert metrics[name] == 0


@pytest.mark.integration
def test_primary_mapped_total_cannot_exceed_raw_reads(tmp_path):
    barcode = "AAAAAAAAAAAAA"
    bam = _write_bam(
        tmp_path / "S1_Aligned.out.bam",
        [(f"read1_{barcode}_", 0, 10)],
    )
    adata = _make_adata(
        [f"S1_{barcode}"],
        [[1]],
        obs_cols={"is_single_cell": [True]},
    )
    h5ad = tmp_path / "S1.raw_feature_bc_matrix.h5ad"
    adata.write_h5ad(h5ad)
    fastp = tmp_path / "S1.R1.preQC.fastp.json"
    fastp.write_text(
        json.dumps({"summary": {"before_filtering": {"total_reads": 0}}}),
        encoding="utf-8",
    )

    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--sample_id",
            "S1",
            "--star_bam",
            str(bam),
            "--raw_count_matrix_h5ad",
            str(h5ad),
            "--fastp_json",
            str(fastp),
            "--barcode_length",
            "13",
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
    )

    assert result.returncode != 0
    assert "Primary alignments with extractable barcodes exceed raw reads" in result.stderr
    assert not (tmp_path / "S1.read_categorization.csv").exists()


@pytest.mark.integration
def test_primary_records_without_extractable_barcodes_cannot_exceed_raw_reads(tmp_path):
    bam = _write_bam(
        tmp_path / "S1_Aligned.out.bam",
        [("read_without_a_barcode", 0, 10)],
    )
    sentinel = tmp_path / "S1.raw_feature_bc_matrix.empty.h5ad"
    sentinel.touch()
    fastp = tmp_path / "S1.R1.preQC.fastp.json"
    fastp.write_text(
        json.dumps({"summary": {"before_filtering": {"total_reads": 0}}}),
        encoding="utf-8",
    )

    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--sample_id",
            "S1",
            "--star_bam",
            str(bam),
            "--raw_count_matrix_h5ad",
            str(sentinel),
            "--fastp_json",
            str(fastp),
            "--barcode_length",
            "13",
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
    )

    assert result.returncode != 0
    assert "Primary BAM records exceed raw reads" in result.stderr
    assert not (tmp_path / "S1.read_categorization.csv").exists()
