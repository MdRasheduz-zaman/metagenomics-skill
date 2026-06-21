import importlib.util
import os
import shutil

import pytest

from metagx import config_builder as cb
from metagx import evidence_pack, tool_advisor

_FIX = os.path.join(os.path.dirname(__file__), "fixtures")
FIXTURE = os.path.join(_FIX, "phylo_tiny.fasta")
DEMO_FIXTURE = os.path.join(_FIX, "phylo_demo.fasta")
_REPO = os.path.dirname(os.path.dirname(__file__))


def _load_phylo_script():
    path = os.path.join(
        os.path.dirname(__file__), "..", "workflow", "scripts", "phylogenetics.py"
    )
    spec = importlib.util.spec_from_file_location("phylogenetics", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _phylo_bin_on_path():
    """Make MAFFT/FastTree/TrimAl resolvable, using the bundled .conda env if present.

    Returns True if the alignment+tree tools are available (so the real-execution
    test can run), else False (CI without the tools — skip cleanly).
    """
    bundled = os.path.join(_REPO, ".conda", "phylogenetics", "bin")
    if os.path.isdir(bundled) and bundled not in os.environ.get("PATH", ""):
        os.environ["PATH"] = bundled + os.pathsep + os.environ.get("PATH", "")
    return bool(shutil.which("mafft") and shutil.which("FastTree"))


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


def test_newick_stats_ignores_internal_support_values():
    """Regression: FastTree writes internal-node support like ')0.993:' — those are
    not leaves. A naive 'token before colon' count wrongly inflated n_leaves."""
    import tempfile
    phylo = _load_phylo_script()
    # 6-taxon tree with internal support values (as FastTree emits)
    tree = ("(taxonB2:0.012,(taxonA3:0.012,(taxonA1:0.0,taxonA2:0.012)0.458:0.0)"
            "0.993:0.067,(taxonB1:0.0,taxonB3:0.012)0.000:0.0);")
    with tempfile.NamedTemporaryFile("w", suffix=".nwk", delete=False) as fh:
        fh.write(tree)
        path = fh.name
    try:
        assert phylo.newick_stats(path)["n_leaves"] == 6
    finally:
        os.unlink(path)


@pytest.mark.skipif(not _phylo_bin_on_path(),
                    reason="MAFFT/FastTree not installed (skips in CI; runs where the env exists)")
def test_phylogenetics_real_execution(tmp_path):
    """End-to-end REAL run of the module: MAFFT -> TrimAl -> FastTree on 6 sequences.

    This is genuine execution (not dry-run): it catches tool-version/parse issues that
    DAG dry-runs cannot. Uses the bundled .conda/phylogenetics env when present.
    """
    phylo = _load_phylo_script()
    out = tmp_path / "phylo"
    out.mkdir()
    aligned = str(out / "aligned.fasta")
    tree = str(out / "tree.nwk")
    jsn = str(out / "phylogenetics.json")
    png = str(out / "tree.png")

    phylo.main(DEMO_FIXTURE, aligned, tree, jsn, png,
               phylo_cfg={"method": "fasttree", "sequence_type": "nt", "trim": True},
               mafft_cfg={"method": "linsi"}, iqtree_cfg={},
               fasttree_cfg={"model": "gtr"}, threads=2)

    import json
    for p in (aligned, tree, jsn, png):
        assert os.path.getsize(p) > 0, f"empty output {p}"
    payload = json.loads(open(jsn).read())
    assert payload["n_sequences"] == 6
    assert payload["tree_method"] == "fasttree"
    assert payload["tree_stats"]["n_leaves"] == 6        # the bug this would have caught
    assert payload["tree_stats"]["total_branch_length"] > 0
    # all six taxa appear in the Newick tree
    newick = open(tree).read()
    for taxon in ("taxonA1", "taxonA2", "taxonA3", "taxonB1", "taxonB2", "taxonB3"):
        assert taxon in newick


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


def test_iqtree_binary_resolves_across_versions(monkeypatch):
    """Version-robustness: resolve iqtree2 / iqtree3 / iqtree, never hardcode one.

    Regression guard for the IQ-TREE 3 break (ships `iqtree3`, not `iqtree2`)."""
    phylo = _load_phylo_script()
    present = {"iqtree3"}
    monkeypatch.setattr(phylo.shutil, "which", lambda c: c if c in present else None)
    assert phylo._iqtree_binary() == "iqtree3"

    present = {"iqtree2", "iqtree3", "iqtree"}   # prefer v2 (the registry's reference)
    assert phylo._iqtree_binary() == "iqtree2"

    present = {"iqtree"}
    assert phylo._iqtree_binary() == "iqtree"

    present = set()                              # none installed -> actionable error
    monkeypatch.setattr(phylo.shutil, "which", lambda c: None)
    with pytest.raises(FileNotFoundError):
        phylo._iqtree_binary()
