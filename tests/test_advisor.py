import json
import os
import tempfile

import pytest

from metagx import advise, catalog, evidence_pack, history, sync_help


def test_recommend_pacbio_clr_has_evidence_sweep():
    rec = evidence_pack.recommend("kraken2", "pacbio_clr")
    assert rec["sweep_suggest"]
    assert max(rec["sweep_suggest"]) <= 0.02
    assert any("avoid" in w.lower() or "validation" in w.lower() for w in rec["warnings"])


def test_recommend_illumina_stable_grid():
    rec = evidence_pack.recommend("kraken2", "illumina")
    assert 0.0 in rec["sweep_suggest"]
    assert rec["sweep_config"]["param"] == "confidence"


def test_registry_warnings_with_high_confidence():
    warns = evidence_pack.registry_warnings(
        "kraken2", "pacbio_clr", {"confidence": 0.1}
    )
    assert warns


def test_parse_help_flags():
    text = """
  -h, --help            Show help
      --confidence FLOAT  Confidence score threshold
"""
    flags = sync_help.parse_help_flags(text)
    names = {f["flag"] for f in flags}
    assert "--help" in names
    assert "--confidence" in names


def test_history_append_and_best(tmp_path):
    hist = tmp_path / "history.jsonl"
    analysis = {
        "platforms": ["illumina"],
        "metrics": {"mean_percent_classified": 80.0},
        "suggestions": [],
        "warnings": [],
        "verdict": "ok",
    }
    cfg = {"project": "t1", "samples": [{"platform": "illumina"}]}
    history.record_from_run(cfg, "c1.yaml", analysis, success=True, path=str(hist))

    analysis2 = dict(analysis)
    analysis2["metrics"] = {"mean_percent_classified": 95.0}
    history.record_from_run(cfg, "c2.yaml", analysis2, success=True, path=str(hist))

    rows = history.read_entries(path=str(hist))
    assert len(rows) == 2
    best = history.best_trial(path=str(hist), metric="mean_percent_classified")
    assert best["metrics"]["mean_percent_classified"] == 95.0


def test_analyze_minimal_config():
    cfg = {
        "project": "empty",
        "outdir": "results",
        "samples": [{"platform": "pacbio_clr", "reads": "x.fq"}],
        "kraken2": {"confidence": 0.1},
    }
    analysis = advise.analyze(cfg)
    assert analysis["verdict"] in ("ok", "marginal", "poor", "good")
    assert "pacbio_clr" in analysis["platforms"]


def test_catalog_lists_evidence():
    cat = catalog.build_catalog()
    assert "kraken2" in cat["registry_tools"]
    assert "kraken2_confidence" in cat["evidence_files"]
    assert cat["planned_modules"] == []
    assert "mafft" in cat["registry_tools"]


def test_write_advisor_outputs(tmp_path):
    cfg = {"project": "adv", "outdir": str(tmp_path / "results")}
    os.makedirs(tmp_path / "results" / "adv", exist_ok=True)
    analysis = {
        "verdict": "ok",
        "warnings": ["test warning"],
        "suggestions": ["try sweep"],
        "metrics": {"mean_percent_classified": 88.0},
        "config_patches": {"sweep": {"param": "confidence", "values": [0.0, 0.01]}},
    }
    paths = advise.write_advisor_outputs(cfg, analysis, advisor_dir=str(tmp_path / "advisor"))
    assert os.path.isfile(paths["advisor_json"])
    with open(paths["advisor_json"]) as fh:
        loaded = json.load(fh)
    assert loaded["verdict"] == "ok"
