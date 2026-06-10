"""Unit tests for the IMRaD paper generator (LaTeX building; no pdflatex needed)."""
import json
import os

from metagx import paper


def test_esc_specials_and_unicode():
    assert paper.esc("a_b") == "a\\_b"
    assert paper.esc("50%") == "50\\%"
    assert paper.esc("A & B") == "A \\& B"
    # unicode that shows up in tool descriptions maps to safe LaTeX
    assert "rightarrow" in paper.esc("C→T")
    assert "geq" in paper.esc("≥5")


def _write(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write(text)


def _cfg(tmp_path, **mods):
    base = {"qc": False, "classify": True, "abundance": True}
    base.update(mods)
    return {
        "project": "p", "outdir": str(tmp_path), "threads": 2,
        "samples": [{"sample": "a1", "r1": "a1.fq", "group": "case"},
                    {"sample": "b1", "r1": "b1.fq", "group": "control"}],
        "db": {"kraken2": "k", "bracken": "k"}, "modules": base,
    }


def test_build_tex_has_imrad_skeleton(tmp_path):
    cfg = _cfg(tmp_path)
    from metagx import report
    manifest = report.build_manifest(cfg)
    tex = paper.build_tex(cfg, manifest)
    assert "\\documentclass" in tex and "\\begin{document}" in tex and "\\end{document}" in tex
    for sec in ("\\section{Introduction}", "\\section{Methods}",
                "\\section{Results}", "\\section{Discussion}", "\\begin{abstract}"):
        assert sec in tex, f"missing {sec}"


def test_generate_writes_tex_and_reads_results(tmp_path):
    cfg = _cfg(tmp_path, stats=True, differential=True)
    outdir = os.path.join(str(tmp_path), "p")
    # a minimal results tree the paper should pick up
    _write(os.path.join(outdir, "summary", "bracken_combined.tsv"),
           "sample\tlabel\tname\tnew_est_reads\tfraction_total_reads\n"
           "a1\tconfidence_0.0\tEscherichia coli\t100\t0.5\n"
           "b1\tconfidence_0.0\tBacteroides_fragilis\t80\t0.4\n")
    _write(os.path.join(outdir, "stats", "differential_abundance.json"),
           json.dumps({"summary": {"group_a": "case", "group_b": "control", "n_a": 1,
                                   "n_b": 1, "n_taxa": 2, "n_permutations": 999, "fdr": 0.05,
                                   "n_significant": 1},
                       "significant_taxa": ["Escherichia coli"]}))
    res = paper.generate(cfg, compile_pdf=False)
    tex_path = res["paths"]["paper_tex"]
    assert os.path.isfile(tex_path)
    tex = open(tex_path).read()
    # taxon names appear, with LaTeX-unsafe underscore escaped
    assert "Escherichia coli" in tex
    assert "Bacteroides\\_fragilis" in tex          # underscore escaped, not raw
    assert "_fragilis" not in tex.replace("\\_", "")  # no raw underscore leaked
    assert "differential" in tex.lower()
    assert res["compiled"] is False                 # we asked not to compile


def test_no_results_is_graceful(tmp_path):
    cfg = _cfg(tmp_path)
    res = paper.generate(cfg, compile_pdf=False)
    tex = open(res["paths"]["paper_tex"]).read()
    assert "run the workflow before generating the paper" in tex.lower()
