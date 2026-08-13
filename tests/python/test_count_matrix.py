"""
Integration test for bin/count_matrix.py.

count_matrix.py builds a sparse cell x gene AnnData from:
  * a barcode_list CSV (cell, io, ioID) - only the io (13bp barcode) column is
    used, as a whitelist of allowed barcodes,
  * a count table / bcGeneSummary (io, gene_id), one row per observed
    (barcode, gene) read,
  * a features tsv (gene_id, gene_name, chromosome).

It groups the count table by (io, gene_id), counts rows, keeps only ios in the
whitelist, and writes out a {sample}.raw_feature_bc_matrix.h5ad. We craft tiny
inputs and assert the resulting per-(barcode, gene) counts and the mito/total
annotations.

NOTE on the barcode_list header: count_matrix.py reads the CSV with the default
header=0 (first row consumed as a header) then renames columns. So the FIRST
row of the barcode_list file is a header row and is NOT a usable barcode. The
fixture below accounts for that with an explicit header line.
"""

import gzip
import os
import subprocess
import sys
import time

import pytest

anndata = pytest.importorskip("anndata")

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
SCRIPT = os.path.join(REPO_ROOT, "bin", "count_matrix.py")


def _run_count_matrix(workdir):
    barcode_list, count_table, features = _write_inputs(workdir)
    return subprocess.run(
        [
            sys.executable, SCRIPT,
            "--barcode_list", str(barcode_list),
            "--count_table", str(count_table),
            "--gene_list", str(features),
            "--sample", "S1",
            "--mixed_species", "False",
            "--mito_chr", "chrM",
        ],
        capture_output=True, text=True, cwd=str(workdir),
    )


def _write_inputs(tmp_path):
    # Whitelist of ios (13bp barcodes). First line is a header (consumed).
    barcode_list = tmp_path / "barcode_list.csv"
    barcode_list.write_text(
        "cell,io,ioID\n"
        "cell1,AAAAAAAAAAAAA,1\n"
        "cell2,CCCCCCCCCCCCC,2\n"
        "cell3,GGGGGGGGGGGGG,3\n"  # whitelisted but never observed
    )

    # bcGeneSummary: one row per observed (io, gene_id) read. Tab separated,
    # no header. Includes:
    #   - barcode AAA..: 2 reads of geneA, 1 read of geneB
    #   - barcode CCC..: 3 reads of geneA
    #   - barcode TTT.. (NOT whitelisted): should be dropped entirely
    count_table = tmp_path / "bcGeneSummary.txt"
    count_table.write_text(
        "AAAAAAAAAAAAA\tENSG_A\n"
        "AAAAAAAAAAAAA\tENSG_A\n"
        "AAAAAAAAAAAAA\tENSG_B\n"
        "CCCCCCCCCCCCC\tENSG_A\n"
        "CCCCCCCCCCCCC\tENSG_A\n"
        "CCCCCCCCCCCCC\tENSG_A\n"
        "TTTTTTTTTTTTT\tENSG_A\n"  # not in whitelist -> dropped
    )

    # Features tsv (gene_id, gene_name, chromosome) WITH header.
    features = tmp_path / "features.tsv"
    features.write_text(
        "gene_id\tgene_name\tchromosome\n"
        "ENSG_A\tGeneA\tchr1\n"
        "ENSG_B\tGeneB\tchrM\n"   # mitochondrial gene
    )
    return barcode_list, count_table, features


@pytest.mark.integration
def test_count_matrix_construction(tmp_path):
    result = _run_count_matrix(tmp_path)
    assert result.returncode == 0, f"script failed:\nSTDOUT:{result.stdout}\nSTDERR:{result.stderr}"

    h5ad = tmp_path / "S1.raw_feature_bc_matrix.h5ad"
    assert h5ad.exists()

    adata = anndata.read_h5ad(str(h5ad))

    # Two whitelisted barcodes were observed (GGG.. never appears in counts).
    obs_names = list(adata.obs_names)
    assert set(obs_names) == {"S1_AAAAAAAAAAAAA", "S1_CCCCCCCCCCCCC"}

    # Genes: GeneA and GeneB present.
    var_names = list(adata.var_names)
    assert set(var_names) == {"GeneA", "GeneB"}

    # Build a {(barcode, gene): count} lookup from the dense matrix.
    dense = adata.X.toarray()
    lut = {}
    for i, bc in enumerate(obs_names):
        for j, gn in enumerate(var_names):
            lut[(bc, gn)] = int(round(dense[i, j]))

    # AAA..: 2 of GeneA, 1 of GeneB.
    assert lut[("S1_AAAAAAAAAAAAA", "GeneA")] == 2
    assert lut[("S1_AAAAAAAAAAAAA", "GeneB")] == 1
    # CCC..: 3 of GeneA, 0 of GeneB.
    assert lut[("S1_CCCCCCCCCCCCC", "GeneA")] == 3
    assert lut[("S1_CCCCCCCCCCCCC", "GeneB")] == 0

    # TTT.. barcode was not whitelisted -> never appears.
    assert not any(bc.endswith("TTTTTTTTTTTTT") for bc in obs_names)

    # Mito annotation: GeneB is on chrM.
    is_mito = dict(zip(adata.var_names, adata.var["is_mito"]))
    assert bool(is_mito["GeneB"]) is True
    assert bool(is_mito["GeneA"]) is False

    # total_counts per barcode: AAA.. -> 3, CCC.. -> 3.
    totals = dict(zip(adata.obs_names, adata.obs["total_counts"]))
    assert totals["S1_AAAAAAAAAAAAA"] == pytest.approx(3.0)
    assert totals["S1_CCCCCCCCCCCCC"] == pytest.approx(3.0)


@pytest.mark.integration
def test_raw_tripartite_archives_are_byte_reproducible(tmp_path):
    """Independent writes at different times produce identical gzip archives."""
    first = tmp_path / "first"
    second = tmp_path / "second"
    first.mkdir()
    second.mkdir()

    first_result = _run_count_matrix(first)
    assert first_result.returncode == 0, first_result.stderr

    # The old gzip.open() implementation stored wall-clock seconds in the
    # matrix header. Crossing a second boundary makes this a regression test
    # for that behavior, in addition to checking the fixed header directly.
    time.sleep(1.1)
    second_result = _run_count_matrix(second)
    assert second_result.returncode == 0, second_result.stderr

    for filename in ("matrix.mtx.gz", "barcodes.tsv.gz", "features.tsv.gz"):
        first_archive = first / filename
        second_archive = second / filename
        assert first_archive.read_bytes() == second_archive.read_bytes(), filename
        with gzip.open(first_archive, "rb") as first_gzip:
            first_payload = first_gzip.read()
        with gzip.open(second_archive, "rb") as second_gzip:
            assert second_gzip.read() == first_payload

    # Bytes 4..7 are the little-endian gzip MTIME field.
    assert (first / "matrix.mtx.gz").read_bytes()[4:8] == b"\x00\x00\x00\x00"


@pytest.mark.integration
def test_count_matrix_empty_count_table(tmp_path):
    """An empty count table -> the script writes an empty .h5ad and exits 0."""
    barcode_list, _, features = _write_inputs(tmp_path)
    empty_counts = tmp_path / "empty.txt"
    empty_counts.write_text("")

    result = subprocess.run(
        [
            sys.executable, SCRIPT,
            "--barcode_list", str(barcode_list),
            "--count_table", str(empty_counts),
            "--gene_list", str(features),
            "--sample", "S1",
            "--mixed_species", "False",
            "--mito_chr", "chrM",
        ],
        capture_output=True, text=True, cwd=str(tmp_path),
    )
    assert result.returncode == 0, f"STDERR:{result.stderr}"
    empty_h5ad = tmp_path / "S1.raw_feature_bc_matrix.empty.h5ad"
    assert empty_h5ad.exists()
    assert empty_h5ad.stat().st_size == 0
