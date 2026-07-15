"""M2: doctor presence-checks the tools THIS config activates, not just the fixed VERSION_FLOORS.

A core-env tool (e.g. blastn for validation) missing is a fail; a --use-conda tool missing is info.
"""
from metagx import doctor


def test_missing_core_config_tool_is_fail(monkeypatch):
    monkeypatch.setattr(doctor.report, "active_tools", lambda cfg: ["blastn"])
    monkeypatch.setattr(doctor.shutil, "which", lambda _exe: None)          # nothing on PATH
    monkeypatch.setattr(doctor.report, "_env_yaml_packages", lambda _p: {"blast": "x"})  # in core
    monkeypatch.setattr(doctor.report, "conda_package", lambda _t: "blast")
    checks = doctor.check_active_tools({"modules": {"validate": True}})
    assert any(c.name == "tool:blastn" and c.status == "fail" for c in checks)


def test_missing_useconda_tool_is_info(monkeypatch):
    monkeypatch.setattr(doctor.report, "active_tools", lambda cfg: ["gtdbtk"])
    monkeypatch.setattr(doctor.shutil, "which", lambda _exe: None)
    monkeypatch.setattr(doctor.report, "_env_yaml_packages", lambda _p: {"blast": "x"})  # gtdbtk NOT core
    monkeypatch.setattr(doctor.report, "conda_package", lambda _t: "gtdbtk")
    checks = doctor.check_active_tools({"modules": {"domain_taxonomy": True}})
    assert any(c.name == "tool:gtdbtk" and c.status == "info" for c in checks)


def test_present_tool_is_ok(monkeypatch):
    monkeypatch.setattr(doctor.report, "active_tools", lambda cfg: ["vsearch"])
    monkeypatch.setattr(doctor.shutil, "which", lambda exe: "/usr/bin/" + exe)
    checks = doctor.check_active_tools({"modules": {}})
    assert any(c.name == "tool:vsearch" and c.status == "ok" for c in checks)


def test_floor_tools_not_duplicated(monkeypatch):
    # kraken2 is in VERSION_FLOORS -> owned by check_tools, not re-reported here
    monkeypatch.setattr(doctor.report, "active_tools", lambda cfg: ["kraken2", "vsearch"])
    monkeypatch.setattr(doctor.shutil, "which", lambda exe: "/usr/bin/" + exe)
    names = {c.name for c in doctor.check_active_tools({"modules": {}})}
    assert "tool:kraken2" not in names and "tool:vsearch" in names


def test_no_cfg_is_noop():
    assert doctor.check_active_tools(None) == []
