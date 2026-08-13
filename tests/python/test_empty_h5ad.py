"""Unit tests for the explicit empty-H5AD sentinel contract."""

from pathlib import Path

import pytest

import empty_h5ad


@pytest.mark.unit
def test_named_readable_zero_byte_h5ad_is_the_only_empty_sentinel(tmp_path):
    sentinel = tmp_path / "sample.raw_feature_bc_matrix.empty.h5ad"
    sentinel.touch()
    ordinary = tmp_path / "sample.raw_feature_bc_matrix.h5ad"
    ordinary.write_bytes(b"nonempty input is parsed by the caller")

    assert empty_h5ad.is_empty_h5ad_sentinel(sentinel) is True
    assert empty_h5ad.is_empty_h5ad_sentinel(ordinary) is False


@pytest.mark.unit
@pytest.mark.parametrize(
    "filename,payload,error",
    [
        ("sample.h5ad", b"", "zero-byte H5AD"),
        (
            "sample.empty.h5ad",
            b"not an empty sentinel",
            "must be zero bytes",
        ),
    ],
)
def test_misnamed_or_nonempty_sentinel_fails(tmp_path, filename, payload, error):
    path = tmp_path / filename
    path.write_bytes(payload)

    with pytest.raises(ValueError, match=error):
        empty_h5ad.is_empty_h5ad_sentinel(path)


@pytest.mark.unit
@pytest.mark.parametrize("kind", ["missing", "directory"])
def test_nonexistent_or_nonregular_h5ad_fails(tmp_path, kind):
    path = tmp_path / "sample.empty.h5ad"
    if kind == "directory":
        path.mkdir()

    with pytest.raises(ValueError, match="existing regular file"):
        empty_h5ad.is_empty_h5ad_sentinel(path)


@pytest.mark.unit
def test_unreadable_h5ad_fails(tmp_path, monkeypatch):
    path = tmp_path / "sample.empty.h5ad"
    path.touch()
    original_open = Path.open

    def deny_target(candidate, *args, **kwargs):
        if candidate == path:
            raise PermissionError("synthetic unreadable input")
        return original_open(candidate, *args, **kwargs)

    monkeypatch.setattr(Path, "open", deny_target)
    with pytest.raises(ValueError, match="not readable"):
        empty_h5ad.is_empty_h5ad_sentinel(path)
