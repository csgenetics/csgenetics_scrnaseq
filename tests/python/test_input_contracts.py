"""Static and behavioral guards for customer input shell boundaries."""

from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


def test_pipeline_centralizes_sample_and_threshold_validation():
    source = (ROOT / "main.nf").read_text(encoding="utf-8")

    assert "def requireSafeSampleId(value)" in source
    assert "[A-Za-z0-9][A-Za-z0-9._-]{0,127}" in source
    assert source.count("requireSafeSampleId(rowList[0])") == 2
    assert "def canonicalManualThreshold(value, label)" in source
    assert "def canonicalMinimumCountThreshold(value)" in source
    assert "params.minimum_count_threshold = canonicalMinimumCountThreshold(" in source
    assert "conflicting manual Cell Caller thresholds" in source


def test_resolved_configuration_uses_utf8_base64_not_raw_shell_json():
    source = (
        ROOT / "modules/local/save_resolved_configuration/main.nf"
    ).read_text(encoding="utf-8")

    assert "json_indented.getBytes('UTF-8').encodeBase64().toString()" in source
    assert "base64 -d > resolved_configuration.txt" in source
    assert "echo '${json_indented}'" not in source
