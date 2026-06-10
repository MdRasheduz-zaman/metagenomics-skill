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


def _tabular_column_counts_consistent(tex: str):
    """Every row in every tabular must have the same cell count as its column spec.

    A row with more '&' than the spec allows is the 'Extra alignment tab' LaTeX fatal
    error — invalid output PDF. This generic check catches header/row drift in ANY table.
    """
    import re
    problems = []
    # spec may contain one level of nested braces (e.g. p{0.45\linewidth}), so allow it
    for m in re.finditer(
            r"\\begin\{tabular\}\{((?:[^{}]|\{[^{}]*\})*)\}(.*?)\\end\{tabular\}", tex, re.S):
        spec, body = m.group(1), m.group(2)
        # count columns in the spec: l/r/c and each p{...} = one column
        n_cols = len(re.findall(r"p\{[^}]*\}|[lrc]", spec))
        for line in body.split("\\\\"):
            line = line.strip()
            if not line or line.startswith("\\hline") or "tabular" in line:
                continue
            cells = line.replace("\\hline", "").count("&") + 1
            if cells > 1 and cells != n_cols:
                problems.append((spec, n_cols, cells, line[:60]))
    return problems


def test_paper_alpha_table_matches_wide_tsv(tmp_path):
    """Regression: alpha_diversity.tsv gained Chao1/ACE/Good's columns (5->8); the paper
    table header was hardcoded to 5, producing invalid LaTeX. Header must track the file."""
    cfg = _cfg(tmp_path, stats=True)
    outdir = os.path.join(str(tmp_path), "p")
    _write(os.path.join(outdir, "summary", "bracken_combined.tsv"),
           "sample\tlabel\tname\tnew_est_reads\tfraction_total_reads\n"
           "a1\tconfidence_0.0\tEscherichia coli\t100\t0.5\n"
           "b1\tconfidence_0.0\tBacteroides fragilis\t80\t0.4\n")
    # the full 8-column alpha table metagx/diversity.py now writes
    _write(os.path.join(outdir, "stats", "alpha_diversity.tsv"),
           "sample\trichness\tchao1\tace\tgoods_coverage\tshannon\tsimpson\tpielou_evenness\n"
           "a1\t29\t29.0\t29.0\t1.0\t3.01\t0.94\t0.89\n"
           "b1\t30\t30.0\t30.6\t0.999\t3.00\t0.94\t0.88\n")
    _write(os.path.join(outdir, "stats", "diversity.json"),
           json.dumps({"n_samples": 2, "n_taxa": 2, "pcoa_explained": [0.9, 0.1],
                       "alpha": [{"sample": "a1"}, {"sample": "b1"}]}))
    res = paper.generate(cfg, compile_pdf=False)
    tex = open(res["paths"]["paper_tex"]).read()
    assert "Chao1" in tex and "Good's cov." in tex          # new columns surfaced
    problems = _tabular_column_counts_consistent(tex)
    assert not problems, f"inconsistent table column counts (invalid LaTeX): {problems}"


def test_no_results_is_graceful(tmp_path):
    cfg = _cfg(tmp_path)
    res = paper.generate(cfg, compile_pdf=False)
    tex = open(res["paths"]["paper_tex"]).read()
    assert "run the workflow before generating the paper" in tex.lower()
