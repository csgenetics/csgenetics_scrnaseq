"""Fail-loud checks for repository scripts executed inside process Conda envs."""

from pathlib import Path

import yaml


REPO_ROOT = Path(__file__).resolve().parents[2]


def _dependencies(environment_name):
    environment = yaml.safe_load(
        (REPO_ROOT / "conda_envs" / environment_name).read_text(encoding="utf-8")
    )
    return environment["dependencies"]


def test_star_conda_environment_supplies_python_for_alignment_count_parser():
    module = (REPO_ROOT / "modules/local/star/main.nf").read_text(encoding="utf-8")
    config = (REPO_ROOT / "conf/conda_envs.config").read_text(encoding="utf-8")
    dependencies = _dependencies("star_samtools.yml")

    assert "star_alignment_counts.py" in module
    assert 'conda_envs/star_samtools.yml' in config
    python_dependencies = [
        dependency for dependency in dependencies
        if isinstance(dependency, str)
        and dependency.strip().split()[0].split("=")[0] == "python"
    ]
    assert python_dependencies
    assert python_dependencies == ["python =3.11"], (
        "STAR's parser runtime must be explicitly pinned, not inherited from a "
        "host or transitive dependency"
    )
