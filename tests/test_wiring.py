"""Wiring-integrity gate: every tool/module must be in sync across all moving parts.

This is the "DAG" guard — it fails when a tool or module is wired into some parts (registry,
config_builder, Snakefile, doctor/dbprovision, advisor, report, MCP, interview docs) but not
others. Adding `metagx/wiring.py` checks; this test enforces zero gaps so a half-wired addition
can't merge green. See wiring.audit() for the individual invariants.
"""
from metagx import config_builder, registry, tool_advisor, wiring


def test_no_wiring_gaps():
    rep = wiring.audit()
    assert rep["ok"], "wiring gaps:\n  - " + "\n  - ".join(rep["gaps"])


def test_audit_detects_an_injected_gap(monkeypatch):
    """The audit must actually catch drift, not vacuously pass — drop a module's advisor entry."""
    patched = dict(tool_advisor.MODULE_TOOLS)
    patched.pop("classify")
    monkeypatch.setattr(tool_advisor, "MODULE_TOOLS", patched)
    rep = wiring.audit()
    assert not rep["ok"]
    assert any("classify" in g and "MODULE_TOOLS" in g for g in rep["gaps"])


def test_every_module_has_advisor_entry():
    for m in config_builder.DEFAULT_MODULES:
        assert m in tool_advisor.MODULE_TOOLS, f"module {m} not in tool_advisor.MODULE_TOOLS"


def test_dbprovision_specs_are_accepted_db_keys():
    for s in dbprovision_specs():
        assert s in config_builder.DB_EXTRA_KEYS, f"db.{s} not accepted by config_builder"


def dbprovision_specs():
    from metagx import dbprovision
    return sorted(dbprovision.SPECS)
