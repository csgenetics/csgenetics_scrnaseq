#!/usr/bin/env python3

"""Parse the alignment-count contract from a STAR ``Log.final.out`` file."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
import re


MAX_SIGNED_64_BIT = (1 << 63) - 1
REQUIRED_FIELDS = {
    "Number of input reads": "input_reads",
    "Uniquely mapped reads number": "uniquely_mapped_reads",
    "Number of reads mapped to multiple loci": "multimapped_reads",
}


@dataclass(frozen=True)
class StarAlignmentCounts:
    """Validated read counts used to route a sample after STAR."""

    input_reads: int
    uniquely_mapped_reads: int
    multimapped_reads: int

    @property
    def aligned_reads(self) -> int:
        return self.uniquely_mapped_reads + self.multimapped_reads


def _parse_count(label: str, raw_value: str) -> int:
    value = raw_value.strip()
    if re.fullmatch(r"[0-9]+", value) is None:
        raise ValueError(
            f"STAR field {label!r} must be a non-negative integer, got {value!r}"
        )

    count = int(value)
    if count > MAX_SIGNED_64_BIT:
        raise ValueError(
            f"STAR field {label!r} exceeds the signed 64-bit limit: {value}"
        )
    return count


def parse_star_final_log(path: str | Path) -> StarAlignmentCounts:
    """Return validated input, unique, multimapper, and total-aligned counts.

    STAR reports unique and accepted multi-locus mappings separately. Both are
    evidence that the sample has alignments and must continue through feature
    assignment. Reads mapped to too many loci are deliberately not included:
    STAR does not emit those alignments under this pipeline's configuration.
    """

    log_path = Path(path)
    values: dict[str, int] = {}

    with log_path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if "|" not in line:
                continue
            raw_label, raw_value = line.split("|", 1)
            label = raw_label.strip()
            attribute = REQUIRED_FIELDS.get(label)
            if attribute is None:
                continue
            if attribute in values:
                raise ValueError(
                    f"STAR field {label!r} occurs more than once in {log_path} "
                    f"(duplicate at line {line_number})"
                )
            values[attribute] = _parse_count(label, raw_value)

    missing = [
        label for label, attribute in REQUIRED_FIELDS.items() if attribute not in values
    ]
    if missing:
        raise ValueError(
            f"STAR log {log_path} is missing required field(s): {', '.join(missing)}"
        )

    counts = StarAlignmentCounts(**values)
    if counts.uniquely_mapped_reads > MAX_SIGNED_64_BIT - counts.multimapped_reads:
        raise ValueError("STAR aligned-read total exceeds the signed 64-bit limit")
    if counts.aligned_reads > counts.input_reads:
        raise ValueError(
            "STAR log is inconsistent: uniquely mapped plus multi-locus reads "
            f"({counts.aligned_reads}) exceeds input reads ({counts.input_reads})"
        )
    return counts


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Extract validated unique and total-aligned read counts from STAR"
    )
    parser.add_argument("log_final_out", help="Path to STAR Log.final.out")
    args = parser.parse_args()

    counts = parse_star_final_log(args.log_final_out)
    print(f"{counts.uniquely_mapped_reads}\t{counts.aligned_reads}")


if __name__ == "__main__":
    main()
