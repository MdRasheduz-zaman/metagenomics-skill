import json
import os
import tempfile

import pytest

from metagx import advise, catalog, evidence_pack, history, sync_help


def test_sample_contexts_from_tsv_sheet_is_shape_complete(tmp_path):
    """Regression: a TSV sample-sheet path must yield contexts with all keys.

    Previously sample_contexts returned a fallback dict missing 'qc_key' for any
    string (path) samples value, crashing the post-run advisor (qc_routing_for) for
    every TSV-based run — the common case.
    """
    from metagx import tool_advisor
    sheet = tmp_path / "samples.tsv"
    sheet.write_text("sample\tr1\tplatform\tlayout\n"
                     "s1\ts1.fastq\tont\tse\n"
                     "s2\ts2.fastq\tillumina\tse\n")
    cfg = {"samples": str(sheet), "modules": {"classify": True}}
    ctxs = tool_advisor.sample_contexts(cfg)
    assert len(ctxs) == 2
    for c in ctxs:                                   # every context is shape-complete
        for key in ("sample", "platform", "library", "layout", "qc_key", "reads"):
            assert key in c
    assert {c["platform"] for c in ctxs} == {"ont", "illumina"}   # real platforms, not fallback
    # the full advisor path must not raise on a TSV-sheet config
    rec = tool_advisor.recommend_config(cfg)
    assert "qc_routing" in rec


def test_diversity_suggestions_low_coverage():
    div = {"n_samples": 3, "core_taxa": [{"taxon": "A"}],
           "alpha": [{"sample": "s1", "goods_coverage": 0.80},
                     {"sample": "s2", "goods_coverage": 0.85},
                     {"sample": "s3", "goods_coverage": 0.99}]}
    sugg, metrics = advise.diversity_suggestions(div)
    assert any("Good's coverage" in s and "deeper sequencing" in s for s in sugg)
    assert metrics["mean_goods_coverage"] < 0.95
    assert metrics["low_coverage_samples"] == ["s1", "s2"]   # the under-sampled ones
    assert metrics["n_core_taxa"] == 1


def test_diversity_suggestions_saturated_with_core_is_quiet():
    div = {"n_samples": 3, "core_taxa": [{"taxon": "A"}, {"taxon": "B"}],
           "alpha": [{"sample": s, "goods_coverage": 0.999} for s in ("s1", "s2", "s3")]}
    sugg, metrics = advise.diversity_suggestions(div)
    assert sugg == []                                         # well-sampled -> no nagging
    assert metrics["n_core_taxa"] == 2


def test_diversity_suggestions_empty_core_flags_heterogeneity():
    div = {"n_samples": 4, "core_taxa": [],
           "alpha": [{"sample": s, "goods_coverage": 0.99} for s in ("a", "b", "c", "d")]}
    sugg, _ = advise.diversity_suggestions(div)
    assert any("Core microbiome is empty" in s for s in sugg)


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
