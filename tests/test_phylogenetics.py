import importlib.util
import os

import pytest

from metagx import config_builder as cb
from metagx import evidence_pack, tool_advisor

FIXTURE = os.path.join(os.path.dirname(__file__), "fixtures", "phylo_tiny.fasta")


def _load_phylo_script():
    path = os.path.join(
        os.path.dirname(__file__), "..", "workflow", "scripts", "phylogenetics.py"
    )
    spec = importlib.util.spec_from_file_location("phylogenetics", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_count_sequences_and_newick_stats():
    import tempfile

    phylo = _load_phylo_script()
    assert phylo.count_sequences(FIXTURE) == 3
    tree = "(s1:0.1,s2:0.2,s3:0.15);"
    with tempfile.NamedTemporaryFile("w", suffix=".nwk", delete=False) as fh:
        fh.write(tree)
        path = fh.name
    try:
        s = phylo.newick_stats(path)
        assert s["n_leaves"] >= 3
        assert s["total_branch_length"] > 0
    finally:
        os.unlink(path)


def test_phylogenetics_config_validation():
    from metagx import registry

    with pytest.raises(registry.ValidationError):
        cb.build_config(
            samples=[{"sample": "s", "r1": "a.fq"}],
            db={"kraken2": "db"},
            modules={"phylogenetics": True},
        )
    cfg = cb.build_config(
        samples=[{"sample": "s", "r1": "a.fq"}],
        db={"kraken2": "db"},
        modules={"phylogenetics": True, "classify": False, "abundance": False},
        phylogenetics={"input": FIXTURE, "method": "iqtree"},
        mafft={"method": "auto"},
        iqtree={"bootstrap": 100},
    )
    assert cfg["modules"]["phylogenetics"] is True
    assert cfg["phylogenetics"]["input"] == FIXTURE
    assert cfg["iqtree"]["bootstrap"] == 100


def test_phylogenetics_evidence_and_optional_module():
    rec = evidence_pack.recommend("mafft", "default", param="method")
    assert rec["value_suggest"] == "auto"
    rec_boot = evidence_pack.recommend("iqtree", "default", param="bootstrap")
    assert rec_boot["value_suggest"] == 1000

    cfg = {
        "samples": [{"sample": "s", "platform": "illumina"}],
        "modules": {"classify": True},
    }
    full = tool_advisor.recommend_config(cfg)
    mods = {o["module"]: o for o in full["optional_modules"]}
    assert "phylogenetics" in mods
    assert mods["phylogenetics"]["ready"] is True


def test_cutadapt_and_metaspades_evidence():
    cut = evidence_pack.recommend("cutadapt", "illumina", param="minimum_length")
    assert cut["value_suggest"] == 100
    mem = evidence_pack.recommend("metaspades", "illumina", param="memory_gb")
    assert mem["value_suggest"] == 250
