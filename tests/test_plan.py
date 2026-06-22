"""Intent-first funnel backend (metagx plan): a goal/preset -> modules + DB checklist.

The whole point is that DB needs fall out of the research-question conversation. These tests
pin the routing: classifier DB surfaces for classifying modules, module DBs come from
dbprovision (single source of truth), each entry tells the LLM how to satisfy it, `have`
marks a DB resolved (not re-asked), and GTDB-Tk carries its never-auto-fetch warning.
"""
import pytest

from metagx import dbprovision, plan


def _names(p):
    return [e["name"] for e in p["databases"]]


def test_classifier_db_surfaces_for_classifying_modules():
    p = plan.plan(modules={"classify": True, "abundance": True})
    kraken = next(e for e in p["databases"] if e["name"] == "kraken2")
    assert kraken["kind"] == "classifier"
    assert kraken["config_key"] == "db.kraken2"
    # all three routes are spelled out for the LLM
    assert "fetch-db" in kraken["fetch"] and "db.build" in kraken["build"]
    assert any("kraken2 database" in q for q in p["questions"])


def test_no_classifier_db_when_not_classifying():
    p = plan.plan(modules={"domain_taxonomy": True}, domains=["viral"])
    assert "kraken2" not in _names(p)


def test_module_dbs_come_from_dbprovision_mapping():
    # plan must not duplicate the module->DB mapping; it delegates to dbprovision.needed_dbs.
    p = plan.plan(modules={"classify": True, "domain_taxonomy": True, "functional": True},
                  domains=["viral"], functional=["annotation", "amr"])
    cfg = {"modules": {"classify": True, "domain_taxonomy": True, "functional": True},
           "domains": ["viral"], "functional": {"annotation": True, "amr": True}}
    expected = set(dbprovision.needed_dbs(cfg)) | {"kraken2"}
    assert set(_names(p)) == expected


def test_preset_seeds_modules_and_summary():
    p = plan.plan(preset="amr-surveillance", functional=["annotation", "amr"])
    assert p["preset"] == "amr-surveillance"
    assert p["summary"]  # when_to_use carried through
    assert "kraken2" in _names(p)  # amr-surveillance classifies
    assert {"bakta", "eggnog", "amrfinderplus"} <= set(_names(p))


def test_have_marks_db_resolved_and_drops_question():
    p = plan.plan(modules={"domain_taxonomy": True}, domains=["prokaryote"], have=["checkm2"])
    by_name = {e["name"]: e for e in p["databases"]}
    assert by_name["checkm2"]["resolved"] is True
    assert by_name["gtdbtk"]["resolved"] is False
    # the resolved DB is not re-asked; the unresolved one is
    assert not any("checkm2" in q for q in p["questions"])
    assert any("gtdbtk" in q for q in p["questions"])


def test_gtdbtk_carries_never_auto_fetch_warning():
    p = plan.plan(modules={"domain_taxonomy": True}, domains=["prokaryote"])
    gtdbtk = next(e for e in p["databases"] if e["name"] == "gtdbtk")
    assert "warn" in gtdbtk and "110 GB" in gtdbtk["warn"]
    assert any("⚠" in q and "gtdbtk" in q for q in p["questions"])


def test_manual_db_routes_to_docs():
    p = plan.plan(modules={"domain_taxonomy": True}, domains=["eukaryote"])
    eukcc = next(e for e in p["databases"] if e["name"] == "eukcc")
    assert eukcc["manual"] is True
    assert "manual download" in eukcc["fetch"]


def test_self_gating_db_is_noted():
    p = plan.plan(modules={"functional": True}, functional=["amr"])
    amr = next(e for e in p["databases"] if e["name"] == "amrfinderplus")
    assert "self-gates" in amr.get("note", "")


def test_unknown_preset_raises_keyerror():
    with pytest.raises(KeyError):
        plan.plan(preset="does-not-exist")
