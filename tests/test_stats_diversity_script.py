"""End-to-end test of the stats/diversity Snakemake script (loaded by path).

Exercises the full main() — matrix build, alpha (incl. chao1/ace/goods), beta
(Bray-Curtis + Jaccard), PCoA, rarefaction, and core microbiome — on a synthetic
combined Bracken table, so the wiring is verified, not just dry-run.
"""
import csv
import importlib.util
import json
from pathlib import Path

_SCRIPTS = Path(__file__).resolve().parents[1] / "workflow" / "scripts"


def _load(name):
    spec = importlib.util.spec_from_file_location(name, _SCRIPTS / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _write_bracken(path):
    rows = [
        # sample, name, new_est_reads, label  — 3 samples, shared + unique taxa
        ("s1", "Bacteroides", 500, "confidence_0.0"),
        ("s1", "Prevotella", 300, "confidence_0.0"),
        ("s1", "Escherichia", 1, "confidence_0.0"),
        ("s2", "Bacteroides", 450, "confidence_0.0"),
        ("s2", "Prevotella", 2, "confidence_0.0"),
        ("s2", "Faecalibacterium", 200, "confidence_0.0"),
        ("s3", "Bacteroides", 600, "confidence_0.0"),
        ("s3", "Faecalibacterium", 150, "confidence_0.0"),
        ("s3", "Akkermansia", 50, "confidence_0.0"),
    ]
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["sample", "name", "new_est_reads", "label"])
        w.writerows(rows)


def test_diversity_script_end_to_end(tmp_path):
    dv_script = _load("diversity")
    bracken = tmp_path / "bracken_combined.tsv"
    _write_bracken(bracken)
    out = tmp_path / "stats"

    paths = {k: str(out / v) for k, v in {
        "alpha": "alpha_diversity.tsv", "beta": "beta_braycurtis.tsv",
        "pcoa": "pcoa.tsv", "json": "diversity.json",
        "barplot": "composition_barplot.png", "pcoa_png": "pcoa.png",
        "rarefaction": "rarefaction.tsv", "rarefaction_png": "rarefaction.png",
        "core": "core_taxa.tsv", "jaccard": "beta_jaccard.tsv",
    }.items()}

    dv_script.main(str(bracken), "confidence_0.0", str(out),
                   paths["alpha"], paths["beta"], paths["pcoa"], paths["json"],
                   paths["barplot"], paths["pcoa_png"], paths["rarefaction"],
                   paths["rarefaction_png"], paths["core"], paths["jaccard"],
                   core_prevalence=0.8)

    # every declared output exists
    for p in paths.values():
        assert Path(p).is_file(), f"missing output {p}"

    # alpha table carries the new richness-estimator columns
    with open(paths["alpha"]) as fh:
        header = next(csv.reader(fh, delimiter="\t"))
    for col in ("richness", "chao1", "ace", "goods_coverage", "shannon"):
        assert col in header

    # core microbiome: Bacteroides is in all 3 samples -> core at 0.8
    with open(paths["core"]) as fh:
        core_taxa = [r["taxon"] for r in csv.DictReader(fh, delimiter="\t")]
    assert "Bacteroides" in core_taxa

    # JSON has the new sections and a rarefaction curve per sample
    data = json.loads(Path(paths["json"]).read_text())
    assert data["n_samples"] == 3
    assert "core_taxa" in data and "rarefaction" in data
    assert set(data["rarefaction"]) == {"s1", "s2", "s3"}
