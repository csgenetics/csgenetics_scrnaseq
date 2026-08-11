"""Release-version consistency checks for customer-facing metadata."""

import re
import shlex
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
SEMVER_PATTERN = r"(?:0|[1-9]\d*)\.(?:0|[1-9]\d*)\.(?:0|[1-9]\d*)"


def release_version() -> str:
    lines = (REPO_ROOT / "version").read_text(encoding="utf-8").splitlines()
    assert len(lines) == 1, "version must contain exactly one line"
    version = lines[0].strip()
    assert re.fullmatch(SEMVER_PATTERN, version), (
        f"version must contain one stable semantic version, found {version!r}"
    )
    return version


def nextflow_manifest_version() -> str:
    config = (REPO_ROOT / "nextflow.config").read_text(encoding="utf-8")
    manifest = re.search(
        r"(?ms)^\s*manifest\s*\{(?P<body>.*?)^\s*\}",
        config,
    )
    assert manifest is not None, "nextflow.config must contain a manifest block"

    assignments = re.findall(
        r"(?m)^\s*version\s*=\s*(['\"])([^'\"]+)\1\s*(?://.*)?$",
        manifest.group("body"),
    )
    assert len(assignments) == 1, (
        "nextflow.config manifest must contain exactly one quoted version assignment"
    )
    return assignments[0][1]


def markdown_fenced_commands(markdown: str):
    """Yield complete shell command lines from fenced Markdown code blocks."""

    in_fence = False
    continued = ""
    for raw_line in markdown.splitlines():
        if raw_line.lstrip().startswith("```"):
            if in_fence:
                assert not continued, "README contains an unfinished shell continuation"
            in_fence = not in_fence
            continue
        if not in_fence:
            continue

        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        if line.endswith("\\"):
            continued += line[:-1].rstrip() + " "
            continue
        yield continued + line
        continued = ""

    assert not in_fence, "README contains an unclosed fenced code block"
    assert not continued, "README contains an unfinished shell continuation"


def readme_remote_launch_versions() -> list[str]:
    readme = (REPO_ROOT / "README.md").read_text(encoding="utf-8")
    versions = []
    target = "csgenetics/csgenetics_scrnaseq"

    for command in markdown_fenced_commands(readme):
        if target not in command:
            continue
        tokens = shlex.split(command, comments=True)
        if len(tokens) < 3 or tokens[:2] != ["nextflow", "run"] or tokens[2] != target:
            continue
        assert "-r" in tokens, (
            f"README remote launch command must pin a release with -r: {command!r}"
        )
        revision_index = tokens.index("-r")
        assert revision_index + 1 < len(tokens), (
            f"README remote launch command has no value after -r: {command!r}"
        )
        versions.append(tokens[revision_index + 1])

    assert versions, (
        "README must contain a fenced remote launch example for "
        "csgenetics/csgenetics_scrnaseq"
    )
    return versions


def changelog_release_version() -> str:
    changelog = (REPO_ROOT / "CHANGELOG.md").read_text(encoding="utf-8")
    release_headings = [
        line.strip() for line in changelog.splitlines() if line.startswith("## ")
    ]
    assert release_headings, "CHANGELOG.md must contain a level-two release heading"

    heading = re.fullmatch(
        rf"## (?P<version>{SEMVER_PATTERN}) - \d{{4}}-\d{{2}}-\d{{2}}",
        release_headings[0],
    )
    assert heading is not None, (
        "the first CHANGELOG.md release heading must be '## X.Y.Z - YYYY-MM-DD', "
        f"found {release_headings[0]!r}"
    )
    return heading.group("version")


def test_version_file_is_a_stable_semantic_version():
    assert release_version()


def test_nextflow_manifest_matches_release_version():
    assert nextflow_manifest_version() == release_version()


def test_readme_remote_launch_examples_match_release_version():
    versions = readme_remote_launch_versions()
    assert set(versions) == {release_version()}, (
        f"README remote launch versions do not agree: {sorted(set(versions))}"
    )


def test_changelog_heading_matches_release_version():
    assert changelog_release_version() == release_version()
