"""
Test for bin/features_names.py (GTF -> features-names tsv).

features_names.py is a monolithic script that reads sys.argv at module import
time, so it cannot be imported. We drive it as a subprocess on a tiny crafted
GTF and assert the resulting tsv exactly.

IMPORTANT - environment: features_names.py relies on gtfparse 1.2.1, which
returns a *pandas* DataFrame, and on the pandas <2.0 string API
(``.str.replace(regex=True)``). Newer, polars-based gtfparse releases break
this script. The main test env carries a modern pandas/anndata stack, so this
single test runs the script in a dedicated pinned env that mirrors
``conda_envs/gtfparse.yml`` (tests/requirements/gtfparse). If that env is not
installed the test is skipped with a clear message.

Behaviour under test (from the script):
  * keeps gene_id, gene_name, seqname (renamed to "chromosome")
  * strips a trailing ".<version>" from gene_id
  * de-duplicates rows
  * replaces NA / empty gene_name with the gene_id
"""

import os
import shutil
import subprocess

import pandas as pd
import pytest

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
SCRIPT = os.path.join(REPO_ROOT, "bin", "features_names.py")
IO_COUNT_EXTRACT = os.path.join(REPO_ROOT, "bin", "io_count_extract")
GTFPARSE_ENV_DIR = os.path.join(REPO_ROOT, "tests", "requirements", "gtfparse")
GTFPARSE_ENV_PYTHON = os.path.join(
    GTFPARSE_ENV_DIR, ".pixi", "envs", "default", "bin", "python"
)


def _gtfparse_python():
    """Return the interpreter for the pinned gtfparse env, or skip."""
    if os.path.exists(GTFPARSE_ENV_PYTHON):
        return [GTFPARSE_ENV_PYTHON]
    if shutil.which("pixi") is not None and os.path.exists(
        os.path.join(GTFPARSE_ENV_DIR, "pixi.toml")
    ):
        return ["pixi", "run", "--manifest-path",
                os.path.join(GTFPARSE_ENV_DIR, "pixi.toml"), "python"]
    pytest.skip(
        "pinned gtfparse test env not installed; run "
        "`pixi install` in tests/requirements/gtfparse"
    )

# A tiny GTF covering:
#  - gene with a version suffix on gene_id (ENSG001.3 -> ENSG001)
#  - gene with a real gene_name
#  - gene with NO gene_name attribute -> gene_name should become gene_id
#  - a duplicate gene line (must be de-duplicated)
TINY_GTF = '\n'.join([
    'chr1\tsrc\tgene\t1\t100\t.\t+\t.\tgene_id "ENSG001.3"; gene_name "Alpha";',
    'chr1\tsrc\ttranscript\t1\t100\t.\t+\t.\tgene_id "ENSG001.3"; gene_name "Alpha";',
    'chr2\tsrc\tgene\t1\t100\t.\t+\t.\tgene_id "ENSG002"; gene_name "Beta";',
    'chr2\tsrc\tgene\t1\t100\t.\t+\t.\tgene_id "ENSG002"; gene_name "Beta";',
    'chrM\tsrc\tgene\t1\t100\t.\t+\t.\tgene_id "ENSG003.1";',
    '',
])


@pytest.mark.integration
def test_features_names_tsv(tmp_path):
    gtf = tmp_path / "tiny.gtf"
    gtf.write_text(TINY_GTF)
    out = tmp_path / "features_names.tsv"

    result = subprocess.run(
        _gtfparse_python() + [SCRIPT, str(gtf), str(out)],
        capture_output=True,
        text=True,
        cwd=str(tmp_path),
    )
    assert result.returncode == 0, f"script failed:\nSTDOUT:{result.stdout}\nSTDERR:{result.stderr}"
    assert out.exists()

    df = pd.read_csv(out, sep="\t")

    # Columns and rename.
    assert list(df.columns) == ["gene_id", "gene_name", "chromosome"]

    rows = {r.gene_id: (r.gene_name, r.chromosome) for r in df.itertuples()}

    # Version suffix stripped, gene_name kept.
    assert "ENSG001" in rows
    assert rows["ENSG001"] == ("Alpha", "chr1")

    # Plain gene with name, de-duplicated to a single row.
    assert rows["ENSG002"] == ("Beta", "chr2")
    assert (df["gene_id"] == "ENSG002").sum() == 1

    # Missing gene_name replaced by gene_id; version suffix stripped.
    assert "ENSG003" in rows
    assert rows["ENSG003"] == ("ENSG003", "chrM")


@pytest.mark.integration
def test_custom_reference_gene_ids_match_io_count_output(tmp_path):
    """GTF feature IDs and XT tags use the same normalization contract."""
    source_ids = [
        "gene-alpha",
        "gene:beta",
        "gene.with.words",
        "gene_under_score",
        "12345",
        "versioned-gene.7",
        "gène-δ",
    ]
    expected_ids = [
        "gene-alpha",
        "gene:beta",
        "gene.with.words",
        "gene_under_score",
        "12345",
        "versioned-gene",
        "gène-δ",
    ]
    gtf = tmp_path / "custom.gtf"
    gtf.write_text("\n".join(
        f'chr1\tsrc\tgene\t{i}\t{i + 1}\t.\t+\t.\t'
        f'gene_id "{gene_id}"; gene_name "Gene{i}";'
        for i, gene_id in enumerate(source_ids, start=1)
    ) + "\n")
    features = tmp_path / "custom_features_names.tsv"

    features_result = subprocess.run(
        _gtfparse_python() + [SCRIPT, str(gtf), str(features)],
        capture_output=True,
        text=True,
        cwd=str(tmp_path),
    )
    assert features_result.returncode == 0, features_result.stderr
    feature_ids = pd.read_csv(features, sep="\t")["gene_id"].tolist()
    assert feature_ids == expected_ids

    sam = "\n".join(
        "read_ACGTACGTACGTA_\t0\tchr1\t1\t255\t1M\t*\t0\t0\tA\tI\t"
        f"XT:Z:{gene_id}"
        for gene_id in source_ids
    ) + "\n"
    extract_result = subprocess.run(
        [IO_COUNT_EXTRACT],
        input=sam,
        capture_output=True,
        text=True,
        cwd=str(tmp_path),
    )
    assert extract_result.returncode == 0, extract_result.stderr
    extracted_ids = [
        line.split("\t", 1)[1]
        for line in extract_result.stdout.splitlines()
    ]
    assert extracted_ids == expected_ids
    assert set(extracted_ids) == set(feature_ids)


@pytest.mark.integration
def test_version_normalization_collision_fails_loudly(tmp_path):
    gtf = tmp_path / "collision.gtf"
    gtf.write_text("\n".join([
        'chr1\tsrc\tgene\t1\t2\t.\t+\t.\tgene_id "custom.1"; gene_name "One";',
        'chr1\tsrc\tgene\t3\t4\t.\t+\t.\tgene_id "custom"; gene_name "Two";',
        '',
    ]))
    output = tmp_path / "features_names.tsv"

    result = subprocess.run(
        _gtfparse_python() + [SCRIPT, str(gtf), str(output)],
        capture_output=True,
        text=True,
        cwd=str(tmp_path),
    )

    assert result.returncode != 0
    assert "collide after version normalization" in result.stderr
    assert "custom.1" in result.stderr
    assert "custom" in result.stderr
    assert not output.exists()
