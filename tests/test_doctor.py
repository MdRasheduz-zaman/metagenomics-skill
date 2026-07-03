"""`metagx doctor` preflight (metagx/doctor.py).

Turns the macOS/arm64 + bioconda landmines into machine-checked diagnostics. These tests
drive the checks with controlled inputs (monkeypatched tool versions / env) so they run the
same everywhere, with no bio tools installed.
"""
import os

import pytest

from metagx import doctor


def _by_name(checks):
    return {c.name: c for c in checks}


def test_workflow_check_ok_in_repo():
    c = doctor.check_workflow()
    assert c.status == "ok"
    assert "Snakefile" in c.message


def test_conda_subdir_leak_warns_on_arm64(monkeypatch):
    monkeypatch.setenv("CONDA_SUBDIR", "osx-64")
    monkeypatch.setattr(doctor.platform, "machine", lambda: "arm64")
    c = doctor.check_conda_subdir_leak()
    assert c.status == "warn"
    assert "base env" in (c.remedy or "")


def test_conda_subdir_clean_when_unset(monkeypatch):
    monkeypatch.delenv("CONDA_SUBDIR", raising=False)
    assert doctor.check_conda_subdir_leak().status == "ok"


def test_apple_silicon_info_emitted(monkeypatch):
    monkeypatch.setattr(doctor.platform, "system", lambda: "Darwin")
    monkeypatch.setattr(doctor.platform, "machine", lambda: "arm64")
    names = {c.name for c in doctor.check_platform()}
    assert "apple-silicon" in names


def test_core_tool_missing_is_fail(monkeypatch):
    # No tools on PATH -> every core tool fails, optional tools are info.
    monkeypatch.setattr(doctor.report, "tool_versions",
                        lambda tools: {t: "not found on PATH" for t in tools})
    checks = _by_name(doctor.check_tools())
    assert checks["tool:kraken2"].status == "fail"
    assert checks["tool:checkv"].status == "info"  # optional module


def test_samtools_downgrade_is_fail_with_abricate_remedy(monkeypatch):
    versions = {t: "9.9" for t in doctor.VERSION_FLOORS}
    versions["samtools"] = "samtools 0.1.19"     # the regression trap
    monkeypatch.setattr(doctor.report, "tool_versions", lambda tools: versions)
    c = _by_name(doctor.check_tools())["tool:samtools"]
    assert c.status == "fail"
    assert "abricate" in (c.remedy or "")


def test_tool_meeting_floor_is_ok(monkeypatch):
    versions = {t: f"{v[0]}.{v[1]}" for t, v in doctor.VERSION_FLOORS.items()}
    monkeypatch.setattr(doctor.report, "tool_versions", lambda tools: versions)
    checks = _by_name(doctor.check_tools())
    assert all(c.status == "ok" for c in checks.values())


def test_database_check_missing_path_fails():
    c = doctor.check_database({"kraken2": "/nonexistent/db/path"})
    assert c.status == "fail"
    assert "fetch-db" in (c.remedy or "") or "build-db" in (c.remedy or "")


def test_database_check_no_config_points_to_fetch():
    c = doctor.check_database(None)
    assert c.status == "info"
    assert "fetch-db" in (c.remedy or "")


def test_database_ok_for_built_db(tmp_path):
    db = tmp_path / "kdb"
    db.mkdir()
    (db / "hash.k2d").write_bytes(b"x" * 1024)
    c = doctor.check_database({"kraken2": str(db)})
    assert c.status == "ok"


def test_param_conflict_mpa_style_with_abundance_fails():
    """use_mpa_style + abundance(Bracken) on is individually-valid but jointly broken (3.3)."""
    cfg = {"modules": {"abundance": True}, "kraken2": {"use_mpa_style": True}}
    checks = doctor.check_param_conflicts(cfg)
    assert any(c.name == "conflict:kraken2.use_mpa_style" and c.status == "fail" for c in checks)


def test_param_conflict_silent_when_module_off():
    """No conflict when the incompatible module is disabled."""
    assert doctor.check_param_conflicts(
        {"modules": {"abundance": False}, "kraken2": {"use_mpa_style": True}}) == []


def test_param_conflict_fires_on_default_on_module(tmp_path):
    """abundance defaults ON; a config that omits `modules:` entirely but sets use_mpa_style
    must still trip the conflict — the effective module map, not the literal block (F1)."""
    cfg = {"kraken2": {"use_mpa_style": True}}  # no `modules:` block at all
    checks = doctor.check_param_conflicts(cfg)
    assert any(c.name == "conflict:kraken2.use_mpa_style" and c.status == "fail" for c in checks)


def test_param_conflict_silent_when_flag_unset():
    """No conflict when the flag itself isn't set, even with the module on."""
    assert doctor.check_param_conflicts(
        {"modules": {"abundance": True}, "kraken2": {"confidence": 0.1}}) == []


def test_version_drift_info_when_installed_differs_from_tested():
    """kraken2 registry records tested_version 2.17.1; a probe reporting a different installed
    version yields an INFO drift check pointing at `metagx refresh`."""
    def probe(tool):
        return {"installed": True, "version": "Kraken version 2.18.0"}
    checks = doctor.check_registry_version_drift(probe=probe)
    assert any(c.name == "version-drift:kraken2" and c.status == "info" for c in checks)


def test_version_drift_silent_when_matching_or_absent():
    assert doctor.check_registry_version_drift(
        probe=lambda t: {"installed": True, "version": "Kraken version 2.17.1"}) == []
    # uninstalled tool is silent (PATH is another check's job)
    assert doctor.check_registry_version_drift(
        probe=lambda t: {"installed": False, "version": None}) == []


def test_format_report_summarizes_failures(monkeypatch):
    checks = [doctor.Check("a", "fail", "broke", remedy="fix it"),
              doctor.Check("b", "ok", "fine")]
    out = doctor.format_report(checks)
    assert "fix it" in out
    assert "1 failure" in out
