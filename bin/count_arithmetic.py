"""Exact, sparse arithmetic for raw count matrices.

Raw count matrices are integer observations even though the public H5AD schema
stores ``X`` as float32 for R/Seurat compatibility.  This module validates that
contract before doing any reduction, converts a copy to canonical CSR/int64,
and never materialises the full matrix as a dense array.
"""

from __future__ import annotations

import numpy as np
from scipy import sparse


INT64_MAX = np.iinfo(np.int64).max
MAX_REDUCTION_CHUNK_VALUES = 1_000_000


def _value_chunks(values: np.ndarray, *, buffer_size=MAX_REDUCTION_CHUNK_VALUES):
    """Yield bounded, vector-sized chunks from any NumPy layout."""
    return np.nditer(
        values,
        flags=["external_loop", "buffered", "zerosize_ok"],
        op_flags=["readonly"],
        order="K",
        buffersize=buffer_size,
    )


def _validate_values(values: np.ndarray, *, context: str) -> None:
    """Validate raw scalar counts without changing their representation."""
    if np.issubdtype(values.dtype, np.complexfloating):
        raise ValueError(f"{context} must contain real counts, not complex values")

    if np.issubdtype(values.dtype, np.floating):
        if values.dtype.itemsize > np.dtype(np.float64).itemsize:
            raise TypeError(
                f"{context} uses unsupported extended floating dtype {values.dtype}"
            )
        # Above this consecutive-integer boundary a floating value no longer
        # proves which raw integer produced it.  Reject it even when that
        # particular floating value happens to be an exactly representable
        # even integer.
        exact_integer_limit = 2 ** (np.finfo(values.dtype).nmant + 1)
        for chunk in _value_chunks(values):
            if not np.all(np.isfinite(chunk)):
                raise ValueError(f"{context} contains non-finite counts")
            if np.any(chunk < 0):
                raise ValueError(f"{context} contains negative counts")
            if not np.all(chunk == np.floor(chunk)):
                raise ValueError(f"{context} contains fractional counts")
            if np.any(chunk > exact_integer_limit):
                raise ValueError(
                    f"{context} contains floating counts above the exact-integer "
                    f"range ({exact_integer_limit}) for {values.dtype}"
                )
        return

    if np.issubdtype(values.dtype, np.integer):
        if values.dtype.itemsize > np.dtype(np.int64).itemsize:
            raise TypeError(
                f"{context} uses unsupported extended integer dtype {values.dtype}"
            )
        for chunk in _value_chunks(values):
            if np.issubdtype(values.dtype, np.signedinteger) and np.any(chunk < 0):
                raise ValueError(f"{context} contains negative counts")
            if np.issubdtype(values.dtype, np.unsignedinteger) and np.any(
                chunk > np.uint64(INT64_MAX)
            ):
                raise OverflowError(
                    f"{context} contains counts outside the int64 range"
                )
        return

    raise TypeError(f"{context} must contain numeric real integer counts")


def _guard_total(values: np.ndarray, *, context: str) -> int:
    """Return the exact Python-integer total, failing before int64 reduction."""
    if values.size == 0:
        return 0

    maximum = int(np.max(values))
    if maximum == 0:
        return 0

    # Each NumPy reduction stays within int64 by construction; Python integers
    # then accumulate the chunks without overflow.  The fixed upper bound also
    # prevents a large copy if a dense, non-contiguous view needs buffering.
    safe_chunk_values = max(1, INT64_MAX // maximum)
    buffer_size = min(MAX_REDUCTION_CHUNK_VALUES, safe_chunk_values)
    total = 0
    for chunk in _value_chunks(values, buffer_size=buffer_size):
        total += int(np.sum(chunk, dtype=np.int64))
        if total > INT64_MAX:
            raise OverflowError(f"{context} total exceeds the int64 range")
    return total


def canonical_count_csr(matrix, *, context: str = "raw count matrix") -> sparse.csr_matrix:
    """Validate ``matrix`` and return a canonical CSR/int64 copy.

    Sparse duplicate entries and explicit zeros are accepted as storage
    details.  They are coalesced/removed on the returned copy.  The source
    matrix is never mutated.  Non-negative inputs allow the exact total to be
    used as a pre-reduction proof that no row, column, or duplicate sum can
    overflow int64.
    """
    if sparse.issparse(matrix):
        if len(matrix.shape) != 2:
            raise ValueError(f"{context} must be two-dimensional")
        values = np.asarray(matrix.data)
        _validate_values(values, context=context)
        expected_total = _guard_total(values, context=context)

        canonical = matrix.astype(np.int64, copy=True).tocsr(copy=False)
        canonical.sum_duplicates()
        canonical.eliminate_zeros()
        canonical.sort_indices()
    else:
        values = np.asarray(matrix)
        if values.ndim != 2:
            raise ValueError(f"{context} must be two-dimensional")
        _validate_values(values, context=context)
        expected_total = _guard_total(values, context=context)
        canonical = sparse.csr_matrix(values.astype(np.int64, copy=True))
        canonical.sort_indices()

    # This is a vector reduction over stored entries, not matrix
    # densification.  The pre-check above makes the int64 sum safe.
    canonical_total = int(canonical.data.sum(dtype=np.int64))
    if canonical_total != expected_total:
        raise ValueError(f"{context} could not be represented exactly as int64")
    return canonical


def row_sums(matrix: sparse.csr_matrix) -> np.ndarray:
    """Return exact int64 row totals from a validated canonical matrix."""
    return np.asarray(matrix.sum(axis=1, dtype=np.int64)).reshape(-1)


def row_nnz(matrix: sparse.csr_matrix) -> np.ndarray:
    """Return the number of detected genes in each row."""
    return np.diff(matrix.indptr).astype(np.int64, copy=False)


def exact_median_floor(values: np.ndarray) -> int:
    """Return the legacy integer median without float conversion.

    Summary metrics historically cast NumPy's median to ``int``.  Counts are
    non-negative, so the even-length half-integer case is rounded down.  The
    split-division expression avoids overflow when the two middle values are
    near int64's upper bound.
    """
    values = np.asarray(values, dtype=np.int64)
    if values.ndim != 1 or values.size == 0:
        raise ValueError("an integer median requires at least one value")
    middle = values.size // 2
    if values.size % 2:
        return int(np.partition(values, middle)[middle])
    middle_values = np.partition(values, (middle - 1, middle))
    lower = int(middle_values[middle - 1])
    upper = int(middle_values[middle])
    return lower // 2 + upper // 2 + (lower % 2 + upper % 2) // 2
