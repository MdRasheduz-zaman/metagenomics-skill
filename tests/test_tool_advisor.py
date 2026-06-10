import yaml

from metagx import evidence_pack, tool_advisor


def test_recommend_config_includes_qc_and_bracken():
    cfg = {
        "project": "t",
        "samples": [{"sample": "s1", "platform": "ont", "r1": "missing.fq"}],
        "modules": {"classify": True, "abundance": True, "qc": True},
        "bracken": {"read_length": 150},
    }
    rec = tool_advisor.recommend_config(cfg)
    assert rec["primary_platform"] == "ont"
    assert rec["qc_routing"]
    assert rec["qc_routing"][0]["pipeline"] == ["porechop_abi", "chopper"]
    assert "bracken" in rec
    assert rec["bracken"]["read_length_by_platform"]["ont"]["suggested_read_length"] == 1000


def test_suggest_bracken_from_median():
    rec = tool_advisor.suggest_bracken_read_length("illumina", median_read_len=248, current=150)
    assert rec["suggested_read_length"] == 250
    assert rec.get("warning")


def test_bracken_read_length_evidence():
    rec = evidence_pack.recommend("bracken", "pacbio_hifi", param="read_length")
    assert rec["value_suggest"] == 1000


def test_kraken2_secondary_in_tool_params():
    cfg = {
        "samples": [{"sample": "s", "platform": "pacbio_clr"}],
        "modules": {"classify": True},
        "kraken2": {"minimum_hit_groups": 3},
    }
    rec = tool_advisor.recommend_config(cfg)
    keys = rec["tool_parameters"].keys()
    assert any(k.startswith("kraken2@") for k in keys)


def test_optional_modules_suggests_reconcile():
    cfg = {
        "samples": [{"sample": "s", "platform": "illumina"}],
        "modules": {"classify": True, "assembly": True, "reconcile": False},
    }
    rec = tool_advisor.recommend_config(cfg)
    mods = {o["module"]: o for o in rec["optional_modules"]}
    assert "reconcile" in mods
    assert mods["reconcile"]["ready"] is True


def test_platform_routing_alternatives():
    route = tool_advisor.qc_routing_for({
        "platform": "illumina", "qc_key": "illumina", "library": "wgs", "layout": "pe",
    })
    alts = {a["tool"] for a in route["alternatives"]}
    assert "trimmomatic" in alts
    assert "fastp" in route["pipeline"]
