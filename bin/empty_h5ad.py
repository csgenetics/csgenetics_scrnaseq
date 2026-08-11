#!/usr/bin/env python3

"""Validate the pipeline's explicit zero-byte H5AD sentinel contract."""

from __future__ import annotations

from pathlib import Path


EMPTY_H5AD_SUFFIX = ".empty.h5ad"


def is_empty_h5ad_sentinel(path: str | Path) -> bool:
    """Return whether *path* is the pipeline's intentional empty H5AD.

    Only an existing, readable, regular, zero-byte file whose basename ends in
    ``.empty.h5ad`` is an empty-data sentinel. Every other path is either a
    normal nonempty H5AD (returned as ``False`` for the caller to parse) or a
    contract violation that fails before a corrupt input can be relabelled as
    intentional emptiness.
    """

    h5ad_path = Path(path)
    if not h5ad_path.is_file():
        raise ValueError(f"H5AD input is not an existing regular file: {h5ad_path}")

    try:
        with h5ad_path.open("rb") as handle:
            first_byte = handle.read(1)
    except OSError as exc:
        raise ValueError(f"H5AD input is not readable: {h5ad_path}") from exc

    is_named_sentinel = h5ad_path.name.endswith(EMPTY_H5AD_SUFFIX)
    is_zero_bytes = first_byte == b""

    if is_zero_bytes:
        if not is_named_sentinel:
            raise ValueError(
                "A zero-byte H5AD is valid only when named "
                f"'*{EMPTY_H5AD_SUFFIX}'; got {h5ad_path.name!r}"
            )
        return True

    if is_named_sentinel:
        raise ValueError(
            f"Named empty H5AD sentinel must be zero bytes: {h5ad_path.name!r}"
        )
    return False
