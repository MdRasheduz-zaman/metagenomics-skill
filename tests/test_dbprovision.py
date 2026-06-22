"""Per-tool module-DB provisioning (genomad/checkv/checkm2/gtdbtk/bakta/amrfinderplus):
the dbprovision registry, the doctor presence-gate, and the config db.provision validation.
Actual downloads are tool-gated and not exercised here (the commands are planned with run=False).
"""
import pytest

from metagx import config_builder as cb
from metagx import dbprovision as p
from metagx import doctor, registry


def test_specs_cover_all_module_dbs():
    assert set(p.SPECS) == {
        "genomad", "checkv", "checkm2", "gtdbtk", "bakta", "amrfinderplus",
        "antismash", "humann_nucleotide", "humann_protein", "eggnog", "metaphlan",
        "eukcc", "emu", "blast"}


def test_manual_tools_have_no_auto_downloader():
    for t in ("eukcc", "emu"):
        assert p.SPECS[t].get("manual")
        res = p.provision(t, "/tmp/whatever-not-real", run=False)
        assert res["ran"] is False and res.get("manual") and "manually" in res["note"]
        assert "manual download" in p.fetch_command(t)


def test_needed_dbs_covers_functional_and_bgc():
    cfg = {"modules": {"functional": True, "bgc": True, "domain_taxonomy": True},
           "domains": ["eukaryote"],
           "functional": {"annotation": True, "pathways": True}}
    need = p.needed_dbs(cfg)
    assert need["eggnog"] == "eggnog" and need["antismash"] == "antismash"
    assert need["humann_nucleotide"] == "humann_nucleotide" and need["eukcc"] == "eukcc"


def test_fetch_command_includes_env_for_gtdbtk():
    assert p.fetch_command("genomad", "/d") == "genomad download-database /d"
    assert "GTDBTK_DATA_PATH=/d" in p.fetch_command("gtdbtk", "/d")
    assert "--type light" in p.fetch_command("bakta", "/d")   # space-friendly default


def test_provision_is_idempotent(tmp_path):
    # a marker present => skipped without running the tool
    (tmp_path / "genomad_db").mkdir()
    (tmp_path / "genomad_db" / "genomad_db.source").write_text("x")
    res = p.provision("genomad", str(tmp_path), run=True)
    assert res["ran"] is False and res.get("skipped") == "already present"


def test_provision_plans_command_without_tool(tmp_path):
    res = p.provision("checkm2", str(tmp_path), run=False)
    assert res["ran"] is False and "checkm2 database --download" in res["command"]


def test_unknown_tool_rejected():
    res = p.provision("nope", "/d", run=False)
    assert res["ok"] is False and "no provisioner" in res["error"]


def test_needed_dbs_maps_modules_to_dbs():
    cfg = {"modules": {"domain_taxonomy": True, "functional": True},
           "domains": ["viral", "prokaryote"], "functional": {"amr": True, "annotation": True}}
    need = p.needed_dbs(cfg)
    assert need == {"genomad": "genomad", "checkv": "checkv", "checkm2": "checkm2",
                    "gtdbtk": "gtdbtk", "bakta": "bakta", "eggnog": "eggnog",
                    "amrfinderplus": "amrfinderplus"}
    # nothing needed when the modules are off
    assert p.needed_dbs({"modules": {}}) == {}


def test_doctor_gate_fails_on_missing_module_db():
    cfg = {"modules": {"domain_taxonomy": True}, "domains": ["viral"], "db": {"kraken2": "k"}}
    checks = doctor.check_module_dbs(cfg)
    by = {c.name: c for c in checks}
    assert by["moduledb:genomad"].status == "fail"
    assert "fetch-db --tool genomad" in by["moduledb:genomad"].remedy


def test_doctor_gate_ok_when_db_present(tmp_path):
    gm = tmp_path / "gm"
    (gm / "genomad_db").mkdir(parents=True)
    (gm / "genomad_db" / "genomad_db.source").write_text("x")
    cv = tmp_path / "cv"
    (cv / "checkv-db-v1.5").mkdir(parents=True)
    (cv / "checkv-db-v1.5" / "genome_db").mkdir()
    (cv / "checkv-db-v1.5" / "genome_db" / "checkv_reps.dmnd").write_text("x")
    cfg = {"modules": {"domain_taxonomy": True}, "domains": ["viral"],
           "db": {"genomad": str(gm), "checkv": str(cv)}}
    statuses = {c.name: c.status for c in doctor.check_module_dbs(cfg)}
    assert statuses["moduledb:genomad"] == "ok" and statuses["moduledb:checkv"] == "ok"


def test_amrfinder_missing_db_is_info_not_fail():
    """The amrfinder rule self-gates on its db, so a missing db is a skip, not a crash."""
    cfg = {"modules": {"functional": True}, "functional": {"amr": True}, "db": {"kraken2": "k"}}
    by = {c.name: c.status for c in doctor.check_module_dbs(cfg)}
    assert by["moduledb:amrfinderplus"] == "info"


def test_config_db_provision_validates():
    samples = [{"sample": "s", "r1": "x.fa", "platform": "illumina", "layout": "se", "contigs": "c.fa"}]
    cfg = cb.build_config(project="p", samples=samples,
                          db={"kraken2": "k", "genomad": "/db/gm", "provision": ["genomad"]},
                          modules={"classify": True, "domain_taxonomy": True}, domains=["viral"])
    assert cfg["db"]["provision"] == ["genomad"]
    with pytest.raises(registry.ValidationError):       # unknown provisioner
        cb.build_config(project="p", samples=samples,
                        db={"kraken2": "k", "provision": ["nope"]}, modules={"classify": True})
    with pytest.raises(registry.ValidationError):       # provision tool without its db path
        cb.build_config(project="p", samples=samples,
                        db={"kraken2": "k", "provision": ["genomad"]}, modules={"classify": True})
