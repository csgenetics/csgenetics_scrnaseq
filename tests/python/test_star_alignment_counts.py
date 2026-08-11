"""Validation tests for STAR alignment evidence used by pipeline routing."""

from pathlib import Path
import subprocess
import sys

import pytest

import star_alignment_counts as sac


SCRIPT = Path(__file__).resolve().parents[2] / "bin" / "star_alignment_counts.py"


def _write_log(path, *, input_reads=20, unique=7, multimapped=3, extra=""):
    path.write_text(
        "                          Started job on | Jan 01 00:00:00\n"
        f"                       Number of input reads | {input_reads}\n"
        f"            Uniquely mapped reads number | {unique}\n"
        f"        Number of reads mapped to multiple loci | {multimapped}\n"
        f"{extra}",
        encoding="utf-8",
    )
    return path


@pytest.mark.unit
def test_parser_retains_unique_and_sums_all_aligned_evidence(tmp_path):
    counts = sac.parse_star_final_log(
        _write_log(tmp_path / "Log.final.out", input_reads=20, unique=7, multimapped=3)
    )

    assert counts.input_reads == 20
    assert counts.uniquely_mapped_reads == 7
    assert counts.multimapped_reads == 3
    assert counts.aligned_reads == 10


@pytest.mark.unit
def test_multimapper_only_sample_has_alignment_evidence(tmp_path):
    counts = sac.parse_star_final_log(
        _write_log(tmp_path / "Log.final.out", input_reads=5, unique=0, multimapped=4)
    )

    assert counts.uniquely_mapped_reads == 0
    assert counts.aligned_reads == 4


@pytest.mark.unit
def test_true_zero_alignment_sample_is_valid(tmp_path):
    counts = sac.parse_star_final_log(
        _write_log(tmp_path / "Log.final.out", input_reads=5, unique=0, multimapped=0)
    )
    assert counts.aligned_reads == 0


@pytest.mark.unit
@pytest.mark.parametrize(
    "missing_label",
    sac.REQUIRED_FIELDS,
)
def test_missing_required_field_fails(tmp_path, missing_label):
    lines = [
        "Number of input reads | 20",
        "Uniquely mapped reads number | 7",
        "Number of reads mapped to multiple loci | 3",
    ]
    path = tmp_path / "Log.final.out"
    path.write_text(
        "\n".join(line for line in lines if not line.startswith(missing_label)) + "\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="missing required field"):
        sac.parse_star_final_log(path)


@pytest.mark.unit
@pytest.mark.parametrize("invalid_value", ["-1", "1.5", "NaN", "", "+2"])
def test_non_integer_or_negative_field_fails(tmp_path, invalid_value):
    path = _write_log(
        tmp_path / "Log.final.out", input_reads=20, unique=invalid_value, multimapped=3
    )
    with pytest.raises(ValueError, match="non-negative integer"):
        sac.parse_star_final_log(path)


@pytest.mark.unit
def test_duplicate_required_field_fails(tmp_path):
    path = _write_log(
        tmp_path / "Log.final.out",
        extra="Uniquely mapped reads number | 7\n",
    )
    with pytest.raises(ValueError, match="occurs more than once"):
        sac.parse_star_final_log(path)


@pytest.mark.unit
def test_individual_field_overflow_fails(tmp_path):
    path = _write_log(
        tmp_path / "Log.final.out",
        input_reads=sac.MAX_SIGNED_64_BIT + 1,
        unique=0,
        multimapped=0,
    )
    with pytest.raises(ValueError, match="signed 64-bit limit"):
        sac.parse_star_final_log(path)


@pytest.mark.unit
def test_aligned_total_overflow_fails(tmp_path):
    path = _write_log(
        tmp_path / "Log.final.out",
        input_reads=sac.MAX_SIGNED_64_BIT,
        unique=sac.MAX_SIGNED_64_BIT,
        multimapped=1,
    )
    with pytest.raises(ValueError, match="aligned-read total exceeds"):
        sac.parse_star_final_log(path)


@pytest.mark.unit
def test_aligned_total_cannot_exceed_input_reads(tmp_path):
    path = _write_log(
        tmp_path / "Log.final.out", input_reads=9, unique=7, multimapped=3
    )
    with pytest.raises(ValueError, match="exceeds input reads"):
        sac.parse_star_final_log(path)


@pytest.mark.integration
def test_cli_emits_unique_and_total_aligned_as_tab_separated_integers(tmp_path):
    path = _write_log(
        tmp_path / "Log.final.out", input_reads=20, unique=0, multimapped=3
    )
    result = subprocess.run(
        [sys.executable, str(SCRIPT), str(path)],
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stderr
    assert result.stdout == "0\t3\n"
