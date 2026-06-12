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
