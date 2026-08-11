"""Exact and memory-bounded count-matrix arithmetic contracts."""

from collections import defaultdict
import tracemalloc
from types import SimpleNamespace

import anndata
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

import count_matrix
from count_arithmetic import (
    INT64_MAX,
    _guard_total,
    canonical_count_csr,
    exact_median_floor,
)
import create_consolidated_report
import create_single_sample_report
import summary_statistics


SS = summary_statistics.SummaryStatistics


def _single_adata(matrix):
    return anndata.AnnData(
        X=matrix,
        obs=pd.DataFrame(
            {"is_single_cell": [True, True, False]},
            index=["cell-1", "cell-2", "noise"],
        ),
        var=pd.DataFrame(
            {"is_mito": [False, True, False, False]},
            index=["nuc-a", "mito", "nuc-b", "nuc-c"],
        ),
    )


def _mixed_adata(matrix, *, hsap_cells=(True, False, False), mmus_cells=(False, True, False)):
    single_cells = np.asarray(hsap_cells) | np.asarray(mmus_cells)
    return anndata.AnnData(
        X=matrix,
        obs=pd.DataFrame(
            {
                "is_single_cell": single_cells,
                "is_hsap_cell": hsap_cells,
                "is_mmus_cell": mmus_cells,
                "is_called_cell": single_cells,
            },
            index=["human-cell", "mouse-cell", "noise"],
        ),
        var=pd.DataFrame(
            {
                "is_hsap": [True, True, False, False, False],
                "is_mmus": [False, False, True, True, False],
                "is_mito_hsap": [False, True, False, False, False],
                "is_mito_mmus": [False, False, False, True, False],
            },
            index=["H-nuc", "H-mito", "M-nuc", "M-mito", "other"],
        ),
    )


def _calculate(adata, *, mixed, reads_pre_qc=120):
    stats = SS.__new__(SS)
    stats.anndata = adata
    stats.mixed = mixed
    stats.metrics_dict = defaultdict(dict)
    stats.metrics_dict["Read QC"]["reads_pre_qc"] = (
        "Number of reads pre-QC",
        reads_pre_qc,
        "test fixture",
    )
    stats._prepare_count_matrix()
    if np.any(stats.single_cell_mask):
        stats.calculate_single_cell_stats()
    else:
        stats.set_single_cell_stats_to_zero()
    return stats


@pytest.mark.unit
@pytest.mark.parametrize(
    "matrix,error,match",
    [
        (np.array([[np.nan]]), ValueError, "non-finite"),
        (np.array([[np.inf]]), ValueError, "non-finite"),
        (np.array([[-1]], dtype=np.int64), ValueError, "negative"),
        (np.array([[1.5]]), ValueError, "fractional"),
        (np.array([[1 + 0j]]), ValueError, "complex"),
        (np.array([[2**24 + 2]], dtype=np.float32), ValueError, "exact-integer"),
        (np.array([[2**53 + 2]], dtype=np.float64), ValueError, "exact-integer"),
        (np.array([["1"]], dtype=object), TypeError, "numeric real"),
    ],
)
def test_raw_count_validation_fails_loud(matrix, error, match):
    with pytest.raises(error, match=match):
        canonical_count_csr(matrix)


@pytest.mark.unit
def test_total_overflow_is_guarded_before_sparse_reduction():
    matrix = sparse.csr_matrix(
        (np.array([INT64_MAX, 1], dtype=np.int64), ([0, 0], [0, 1])),
        shape=(1, 2),
    )
    with pytest.raises(OverflowError, match="total exceeds"):
        canonical_count_csr(matrix)


@pytest.mark.unit
def test_integer_median_is_exact_near_int64_limit():
    assert exact_median_floor(np.array([INT64_MAX - 1, INT64_MAX])) == INT64_MAX - 1


@pytest.mark.unit
def test_total_guard_uses_vector_chunks_not_python_scalar_iteration():
    class NoFlatIteration(np.ndarray):
        @property
        def flat(self):
            raise AssertionError("scalar flat iteration is forbidden")

    values = np.ones(2_000_000, dtype=np.int64).view(NoFlatIteration)
    assert _guard_total(values, context="test counts") == 2_000_000


@pytest.mark.unit
def test_uint64_boundary_is_supported_but_out_of_range_value_fails():
    accepted = canonical_count_csr(np.array([[INT64_MAX]], dtype=np.uint64))
    assert accepted.data.tolist() == [INT64_MAX]

    with pytest.raises(OverflowError, match="outside the int64 range"):
        canonical_count_csr(np.array([[INT64_MAX + 1]], dtype=np.uint64))


@pytest.mark.unit
def test_extended_float_dtype_is_rejected_explicitly_when_available():
    if not hasattr(np, "float128"):
        pytest.skip("NumPy has no extended float dtype on this platform")
    with pytest.raises(TypeError, match="unsupported extended floating dtype"):
        canonical_count_csr(np.array([[1]], dtype=np.float128))


@pytest.mark.unit
def test_sparse_canonicalization_coalesces_copy_without_source_mutation():
    source = sparse.csr_matrix(
        (
            np.array([2.0, 0.0, 3.0, 4.0], dtype=np.float32),
            np.array([1, 0, 1, 2], dtype=np.int32),
            np.array([0, 4], dtype=np.int32),
        ),
        shape=(1, 3),
    )
    before = (source.data.copy(), source.indices.copy(), source.indptr.copy())

    canonical = canonical_count_csr(source)

    np.testing.assert_array_equal(source.data, before[0])
    np.testing.assert_array_equal(source.indices, before[1])
    np.testing.assert_array_equal(source.indptr, before[2])
    assert source.nnz == 4
    assert canonical.dtype == np.int64
    assert canonical.has_canonical_format
    assert canonical.nnz == 2
    np.testing.assert_array_equal(canonical.data, [5, 4])
    np.testing.assert_array_equal(canonical.indices, [1, 2])


@pytest.mark.unit
def test_single_species_statistics_are_hand_calculated_exactly():
    matrix = np.array(
        [
            [2, 1, 0, 4],
            [0, 0, 3, 0],
            [9, 9, 9, 9],
        ],
        dtype=np.float32,
    )
    stats = _calculate(_single_adata(matrix), mixed=False)

    assert stats.num_cells == 2
    assert stats.raw_reads_per_cell == 60
    assert stats.mean_total_counts_per_cell == 5
    assert stats.median_total_counts_per_cell == 5
    assert stats.mean_genes_detected_per_cell == 2
    assert stats.median_genes_detected_per_cell == 2
    assert stats.mean_nuclear_genes_detected_per_cell == 1.5
    assert stats.median_nuclear_genes_detected_per_cell == 1
    assert stats.mean_mito_genes_detected_per_cell == 0.5
    assert stats.median_mito_genes_detected_per_cell == 0
    assert stats.percentage_counts_from_mito == 10
    assert stats.num_unique_genes_detected_across_sample == 4
    assert stats.total_genes_detected_across_sample == 4


@pytest.mark.unit
def test_mixed_statistics_preserve_species_matched_and_all_gene_definitions():
    matrix = sparse.csr_matrix(
        np.array(
            [
                [4, 1, 7, 0, 2],
                [6, 0, 3, 2, 1],
                [9, 9, 9, 9, 9],
            ],
            dtype=np.float32,
        )
    )
    stats = _calculate(_mixed_adata(matrix), mixed=True)

    assert stats.num_cells_total == 2
    assert stats.num_cells_Hsap == 1
    assert stats.num_cells_Mmus == 1
    assert stats.mean_total_counts_per_cell_total == 5
    assert stats.mean_total_counts_per_cell_Hsap == 5
    assert stats.mean_total_counts_per_cell_Mmus == 5
    assert stats.mean_genes_detected_per_cell_total == 2
    assert stats.mean_nuclear_genes_detected_per_cell_total == 1
    assert stats.mean_mito_genes_detected_per_cell_total == 1
    assert stats.percentage_counts_from_mito_total == 30
    assert stats.percentage_counts_from_mito_Hsap == 20
    assert stats.percentage_counts_from_mito_Mmus == 40

    # Total-detection metrics intentionally consider every gene in every
    # single cell, including off-species and unclassified genes.
    assert stats.num_unique_genes_detected_across_sample_total == 5
    assert stats.total_genes_detected_across_sample_total == 8
    assert stats.num_unique_genes_detected_across_sample_Hsap == 2
    assert stats.num_unique_genes_detected_across_sample_Mmus == 2


@pytest.mark.unit
def test_missing_species_and_zero_mito_counts_produce_numeric_zeros():
    matrix = sparse.csr_matrix(
        np.array(
            [
                [4, 0, 0, 0, 0],
                [0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0],
            ],
            dtype=np.float32,
        )
    )
    adata = _mixed_adata(
        matrix,
        hsap_cells=(True, False, False),
        mmus_cells=(False, False, False),
    )
    stats = _calculate(adata, mixed=True)

    assert stats.num_cells_total == 1
    assert stats.num_cells_Hsap == 1
    assert stats.num_cells_Mmus == 0
    assert stats.percentage_counts_from_mito_total == 0
    assert stats.percentage_counts_from_mito_Hsap == 0
    for name in (
        "raw_reads_per_cell_Mmus",
        "mean_total_counts_per_cell_Mmus",
        "median_total_counts_per_cell_Mmus",
        "percentage_counts_from_mito_Mmus",
        "num_unique_genes_detected_across_sample_Mmus",
    ):
        assert getattr(stats, name) == 0


@pytest.mark.unit
def test_no_cells_keep_the_complete_zero_metric_contract():
    adata = _single_adata(sparse.csr_matrix((3, 4), dtype=np.float32))
    adata.obs["is_single_cell"] = False
    stats = _calculate(adata, mixed=False)

    for name in (
        "num_cells",
        "raw_reads_per_cell",
        "mean_total_counts_per_cell",
        "median_total_counts_per_cell",
        "mean_genes_detected_per_cell",
        "percentage_counts_from_mito",
        "num_unique_genes_detected_across_sample",
        "total_genes_detected_across_sample",
    ):
        assert getattr(stats, name) == 0


@pytest.mark.unit
@pytest.mark.parametrize(
    "hsap_cells,mmus_cells",
    [
        ((True, False, False), (True, True, False)),
        ((False, False, False), (False, True, False)),
    ],
    ids=["cell-in-both-species", "single-cell-in-neither-species"],
)
def test_inconsistent_mixed_cell_classification_fails(hsap_cells, mmus_cells):
    adata = _mixed_adata(
        sparse.csr_matrix((3, 5), dtype=np.float32),
        hsap_cells=hsap_cells,
        mmus_cells=mmus_cells,
    )
    if not any(hsap_cells) and any(mmus_cells):
        # Make an explicit single cell which is absent from both species masks.
        adata.obs["is_single_cell"] = [True, True, False]
    with pytest.raises(ValueError, match="partitioned exactly once"):
        _calculate(adata, mixed=True)


@pytest.mark.unit
def test_dense_and_sparse_inputs_emit_identical_metric_csv(tmp_path, monkeypatch):
    matrix = np.array(
        [
            [2, 1, 0, 4],
            [0, 0, 3, 0],
            [9, 9, 9, 9],
        ],
        dtype=np.float32,
    )
    outputs = []
    for label, representation in (
        ("dense", matrix),
        ("sparse", sparse.csr_matrix(matrix)),
    ):
        stats = _calculate(_single_adata(representation), mixed=False)
        stats.populate_cell_stats_in_metrics_dict_single_species()
        stats.sample_id = "sample"
        output_dir = tmp_path / label
        output_dir.mkdir()
        monkeypatch.chdir(output_dir)
        stats.write_out_dict_to_csv()
        outputs.append((output_dir / "sample.metrics.csv").read_bytes())

    assert outputs[0] == outputs[1]


@pytest.mark.integration
def test_h5ad_round_trip_uses_exact_sparse_cell_statistics(tmp_path):
    matrix = sparse.csr_matrix(
        np.array(
            [
                [2, 1, 0, 4],
                [0, 0, 3, 0],
                [9, 9, 9, 9],
            ],
            dtype=np.float32,
        )
    )
    h5ad = tmp_path / "sample.raw_feature_bc_matrix.h5ad"
    _single_adata(matrix).write_h5ad(h5ad)

    stats = SS.__new__(SS)
    stats.args = SimpleNamespace(h5ad=str(h5ad))
    stats.mixed = False
    stats.metrics_dict = defaultdict(dict)
    stats.metrics_dict["Read QC"]["reads_pre_qc"] = (
        "Number of reads pre-QC",
        120,
        "test fixture",
    )
    stats.get_cell_stats()

    assert stats.count_matrix.dtype == np.int64
    assert stats.metrics_dict["Cell metrics"]["mean_total_counts_per_cell"][1] == 5
    assert stats.metrics_dict["Cell metrics"]["percentage_counts_from_mito"][1] == 10


@pytest.mark.unit
def test_exact_csv_spelling_keeps_the_two_decimal_report_format():
    old_float32_mean = np.mean(np.array([10, 10, 11], dtype=np.float32))
    exact_mean = 31 / 3
    assert str(old_float32_mean) != str(exact_mean)

    values = (str(old_float32_mean), str(exact_mean))
    assert {
        create_consolidated_report.format_number_to_string(value)
        for value in values
    } == {"10.33"}
    assert {
        create_single_sample_report.SingleSampleHTMLReport.format_number_to_string(value)
        for value in values
    } == {"10.33"}


@pytest.mark.unit
def test_exact_arithmetic_can_correct_the_last_displayed_decimals():
    # Every per-cell total is exactly representable in float32, but the legacy
    # float32 mean loses the half-count during accumulation. The report keeps
    # its two-decimal presentation while truthfully displaying the correction.
    per_cell_totals = np.arange(10, dtype=np.int64) + 10_000_000
    old_float32_mean = np.mean(per_cell_totals.astype(np.float32))
    exact_mean = int(per_cell_totals.sum()) / per_cell_totals.size

    assert create_consolidated_report.format_number_to_string(
        str(old_float32_mean)
    ) == "10000004.00"
    assert create_consolidated_report.format_number_to_string(
        str(exact_mean)
    ) == "10000004.50"
    assert create_single_sample_report.SingleSampleHTMLReport.format_number_to_string(
        str(exact_mean)
    ) == "10000004.50"


@pytest.mark.integration
def test_count_matrix_annotations_use_sparse_original_integer_counts(monkeypatch):
    original_counts = sparse.csr_matrix(
        np.array([[2, 1], [0, 3]], dtype=np.int64)
    )
    # Deliberately make X disagree: annotations must be derived from the
    # original integral matrix, not its compatibility representation.
    adata = anndata.AnnData(
        X=sparse.csr_matrix(np.array([[99, 99], [99, 99]], dtype=np.float32)),
        var=pd.DataFrame(
            {"gene_name": ["a", "b"], "chromosome": ["chr1", "chrM"]},
            index=["a", "b"],
        ),
    )
    instance = count_matrix.CountMatrix.__new__(count_matrix.CountMatrix)
    instance.sparse_matrix = original_counts
    instance.anndata_obj = adata
    instance.mixed_species = False
    instance.mito_chr = "chrM"

    def forbidden(*args, **kwargs):
        raise AssertionError("full-matrix densification is forbidden")

    monkeypatch.setattr(sparse.csr_matrix, "toarray", forbidden)
    monkeypatch.setattr(sparse.csr_matrix, "todense", forbidden)
    instance.annotate_count_matrix()

    np.testing.assert_array_equal(instance.anndata_obj.obs["total_counts"], [3, 3])
    assert instance.anndata_obj.obs["total_counts"].dtype == np.float32


@pytest.mark.unit
def test_count_matrix_rejects_row_totals_outside_float32_obs_contract():
    original_counts = sparse.csr_matrix(
        np.array([[2**24 + 1]], dtype=np.int64)
    )
    instance = count_matrix.CountMatrix.__new__(count_matrix.CountMatrix)
    instance.sparse_matrix = original_counts
    instance.anndata_obj = anndata.AnnData(
        X=original_counts.astype(np.float32),
        var=pd.DataFrame(
            {"gene_name": ["a"], "chromosome": ["chr1"]}, index=["a"]
        ),
    )
    instance.mixed_species = False
    instance.mito_chr = "chrM"

    with pytest.raises(OverflowError, match="exact float32"):
        instance.annotate_count_matrix()


@pytest.mark.integration
def test_large_logical_matrix_stays_sparse_and_memory_bounded(monkeypatch):
    n_obs = 50_000
    n_vars = 30_000
    n_values = 4_000
    positions = np.arange(n_values, dtype=np.int64)
    matrix = sparse.csr_matrix(
        (
            np.ones(n_values, dtype=np.float32),
            (positions % n_obs, (positions * 7_919) % n_vars),
        ),
        shape=(n_obs, n_vars),
    )
    adata = anndata.AnnData(
        X=matrix,
        obs=pd.DataFrame(
            {"is_single_cell": np.ones(n_obs, dtype=bool)},
            index=[f"cell-{i}" for i in range(n_obs)],
        ),
        var=pd.DataFrame(
            {"is_mito": np.arange(n_vars) < 10},
            index=[f"gene-{i}" for i in range(n_vars)],
        ),
    )

    def forbidden(*args, **kwargs):
        raise AssertionError("full-matrix densification is forbidden")

    monkeypatch.setattr(sparse.csr_matrix, "toarray", forbidden)
    monkeypatch.setattr(sparse.csr_matrix, "todense", forbidden)
    monkeypatch.setattr(anndata.AnnData, "to_df", forbidden)

    tracemalloc.start()
    stats = _calculate(adata, mixed=False)
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    assert stats.num_cells == n_obs
    assert stats.total_genes_detected_across_sample == n_values
    assert peak < 64 * 1024 * 1024
