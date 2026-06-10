"""Full IMRaD manuscript generation (LaTeX -> PDF via pdflatex).

The interview already captures the *experimental design* (samples, platforms, libraries,
groups), the *methods* (every tool + exact parameters, from the registries) and, after a run,
the *results* (classification, diversity, differential abundance, reconciliation, ...). This
module elaborates those into a structured Introduction / Methods / Results / Discussion paper
and compiles it to PDF with pdflatex — a publishable first draft, not a fabricated one: every
number is read back from the result files, and interpretation is framed as caveat-aware
discussion for the author to refine.

Reuses ``metagx.report`` for the manifest (tool versions, commands), the Methods paragraph,
and the citation set, so the paper stays in lockstep with the registries.
"""

from __future__ import annotations

import csv
import datetime as _dt
import glob
import json
import os
import shutil
import subprocess
from typing import Any, Dict, List, Optional

from . import __version__, registry, report

# --------------------------------------------------------------------------- #
# LaTeX helpers                                                               #
# --------------------------------------------------------------------------- #
_LATEX_SPECIALS = {
    "\\": r"\textbackslash{}", "&": r"\&", "%": r"\%", "$": r"\$", "#": r"\#",
    "_": r"\_", "{": r"\{", "}": r"\}", "~": r"\textasciitilde{}",
    "^": r"\textasciicircum{}",
}

# Unicode that shows up in tool descriptions / taxon names -> safe LaTeX (default fonts).
_UNICODE = {
    "→": r"$\rightarrow$", "←": r"$\leftarrow$", "≥": r"$\geq$", "≤": r"$\leq$",
    "×": r"$\times$", "±": r"$\pm$", "π": r"$\pi$", "α": r"$\alpha$", "β": r"$\beta$",
    "μ": r"$\mu$", "–": "--", "—": "---", "’": "'", "‘": "'", "“": "``", "”": "''",
    "≈": r"$\approx$", "°": r"$^{\circ}$",
}


def esc(text: Any) -> str:
    """Escape a string for LaTeX body text (handles specials + common unicode)."""
    out = []
    for ch in str(text):
        if ch in _UNICODE:
            out.append(_UNICODE[ch])
        else:
            out.append(_LATEX_SPECIALS.get(ch, ch))
    return "".join(out)


def _fig(path: str, outdir: str, caption: str, label: str, width: str = "0.8") -> str:
    """A figure float referencing an image by path relative to the report dir."""
    rel = os.path.relpath(path, outdir)
    return (
        "\\begin{figure}[H]\n\\centering\n"
        f"\\includegraphics[width={width}\\linewidth]{{{rel}}}\n"
        f"\\caption{{{esc(caption)}}}\n\\label{{fig:{label}}}\n\\end{{figure}}\n"
    )


def _table(header: List[str], rows: List[List[Any]], caption: str, label: str,
           wrap_first: bool = False) -> str:
    # wrap_first: make the first column a wrapping paragraph (for long taxon names) so wide
    # tables don't overflow the text block.
    first = "p{0.45\\linewidth}" if wrap_first else "l"
    cols = first + "r" * (len(header) - 1)
    head = " & ".join(f"\\textbf{{{esc(h)}}}" for h in header) + " \\\\\n\\hline\n"
    body = "".join(" & ".join(esc(c) for c in r) + " \\\\\n" for r in rows)
    return (
        "\\begin{table}[H]\n\\centering\n\\small\n"
        f"\\caption{{{esc(caption)}}}\n\\label{{tab:{label}}}\n"
        f"\\begin{{tabular}}{{{cols}}}\n\\hline\n{head}{body}\\hline\n"
        "\\end{tabular}\n\\end{table}\n"
    )


# --------------------------------------------------------------------------- #
# Design / results readers                                                    #
# --------------------------------------------------------------------------- #
def _design(cfg: Dict[str, Any]) -> Dict[str, Any]:
    recs = report._records(cfg)
    plats = sorted({str(r.get("platform", "illumina")).lower() for r in recs}) or ["illumina"]
    libs = sorted({str(r.get("library", "wgs")).lower() for r in recs}) or ["wgs"]
    groups: Dict[str, int] = {}
    for r in recs:
        g = r.get("group")
        if g:
            groups[str(g)] = groups.get(str(g), 0) + 1
    return {"n_samples": len(recs) or 1, "platforms": plats, "libraries": libs,
            "groups": groups}


def _top_taxa(outdir: str, top_n: int = 12) -> List[Dict[str, str]]:
    path = os.path.join(outdir, "summary", "bracken_combined.tsv")
    if not os.path.isfile(path):
        return []
    with open(path) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    rows.sort(key=lambda r: float(r.get("fraction_total_reads", 0) or 0), reverse=True)
    return rows[:top_n]


def _read_json(path: str) -> Optional[dict]:
    if not os.path.isfile(path):
        return None
    try:
        with open(path) as fh:
            return json.load(fh)
    except (OSError, json.JSONDecodeError):
        return None


# Friendly column labels for the alpha-diversity table (any unlisted column is
# title-cased). Derived from the TSV header so the table can never drift out of sync
# with the columns metagx/diversity.py actually writes.
_ALPHA_LABELS = {
    "sample": "Sample", "richness": "Richness", "chao1": "Chao1", "ace": "ACE",
    "goods_coverage": "Good's cov.", "shannon": "Shannon", "simpson": "Simpson",
    "pielou_evenness": "Evenness",
}


def _alpha_table(outdir: str) -> "Optional[tuple[List[str], List[List[str]]]]":
    """(header, rows) for the alpha-diversity table, read straight from the TSV.

    The header tracks the file's columns, so adding a metric to alpha_diversity.tsv
    flows into the paper with no change here (previously a hardcoded 5-column header
    desynced from the 8-column file and produced invalid LaTeX)."""
    path = os.path.join(outdir, "stats", "alpha_diversity.tsv")
    if not os.path.isfile(path):
        return None
    with open(path) as fh:
        reader = csv.reader(fh, delimiter="\t")
        cols = next(reader, [])
        if not cols:
            return None
        header = [_ALPHA_LABELS.get(c, c.replace("_", " ").title()) for c in cols]
        rows = [r for r in reader if r]
    return header, rows


# --------------------------------------------------------------------------- #
# Sections                                                                    #
# --------------------------------------------------------------------------- #
def _title(cfg: Dict[str, Any], design: Dict[str, Any]) -> str:
    mods = cfg.get("modules", {})
    aim = "Taxonomic profiling"
    if mods.get("differential"):
        aim = "Comparative taxonomic profiling"
    elif mods.get("assembly"):
        aim = "Assembly-based metagenomic characterization"
    libs = "/".join(design["libraries"])
    plats = ", ".join(design["platforms"])
    return f"{aim} of the {esc(cfg.get('project', 'study'))} {libs} dataset ({plats})"


def _abstract(cfg: Dict[str, Any], design: Dict[str, Any], outdir: str) -> str:
    mods = cfg.get("modules", {})
    s = [f"We analysed {design['n_samples']} sequencing sample(s) "
         f"({', '.join(design['platforms'])}; {', '.join(design['libraries'])} library) "
         f"with metagx, a registry-driven metagenomics workflow."]
    cm = report.classification_metrics(outdir)
    if cm:
        vals = [m.get("percent_classified") for m in cm.values() if m.get("percent_classified") is not None]
        if vals:
            s.append(f"Across reports, {min(vals):.1f}--{max(vals):.1f}\\% of reads were "
                     f"taxonomically classified.")
    da = _read_json(os.path.join(outdir, "stats", "differential_abundance.json"))
    if da:
        sm = da.get("summary", {})
        s.append(f"Differential-abundance testing between {esc(sm.get('group_a','A'))} and "
                 f"{esc(sm.get('group_b','B'))} identified {sm.get('n_significant',0)} taxa at "
                 f"FDR {sm.get('fdr','')}.")
    if mods.get("assembly"):
        s.append("De novo assembly and downstream characterization were also performed.")
    s.append("All parameters, tool versions and exact commands are reported for reproducibility.")
    return " ".join(s)


def _introduction(cfg: Dict[str, Any], design: Dict[str, Any]) -> str:
    mods = cfg.get("modules", {})
    p = [
        "Metagenomic sequencing characterizes microbial communities directly from a sample, "
        "without cultivation, by classifying reads against reference databases and, where "
        "depth allows, reconstructing genomes de novo. Robust analysis requires careful "
        "quality control, an appropriately complete reference, and statistics that respect the "
        "compositional nature of sequencing data."
    ]
    aim = ["The aim of this study was to"]
    goals = []
    if mods.get("classify"):
        goals.append("profile the taxonomic composition of the samples")
    if cfg.get("sweep"):
        goals.append("assess the sensitivity of classification to the confidence threshold")
    if mods.get("differential"):
        gs = design["groups"]
        goals.append(f"identify taxa that differ between sample groups ({', '.join(esc(g) for g in gs)})")
    if mods.get("assembly"):
        goals.append("assemble and characterize the recovered genomic content")
    if mods.get("functional") or mods.get("bgc"):
        goals.append("survey functional and biosynthetic potential")
    if not goals:
        goals = ["characterize the sequenced material"]
    aim.append(", ".join(goals) + ".")
    p.append(" ".join(aim))
    return "\n\n".join(p)


def _methods(cfg: Dict[str, Any], manifest: Dict[str, Any], design: Dict[str, Any]) -> str:
    p = []
    g = ", ".join(f"{esc(k)} (n={v})" for k, v in design["groups"].items())
    design_line = (
        f"The dataset comprised {design['n_samples']} sample(s) sequenced on "
        f"{', '.join(design['platforms'])} ({', '.join(design['libraries'])} library)."
    )
    if g:
        design_line += f" Samples were assigned to groups: {g}."
    p.append(design_line)
    # the registry-derived Methods paragraph (tools + parameters, version-aware).
    # It is plain prose (built for markdown) -> escape it before placing in LaTeX.
    p.append(esc(report.methods_paragraph(cfg, manifest)))
    return "\n\n".join(p)


def _results(cfg: Dict[str, Any], manifest: Dict[str, Any], outdir: str) -> str:
    mods = cfg.get("modules", {})
    repdir = os.path.join(outdir, "report")   # paper.tex lives here; figures are relative to it
    blocks: List[str] = []

    cm = manifest.get("metrics", {}).get("classification", {})
    if cm:
        rows = [[k, m.get("percent_classified", "")] for k, m in cm.items()]
        blocks.append("\\subsection*{Classification}")
        vals = [m.get("percent_classified") for m in cm.values() if m.get("percent_classified") is not None]
        if vals:
            blocks.append(
                f"The fraction of classified reads ranged from {min(vals):.1f}\\% to "
                f"{max(vals):.1f}\\% across the reports analysed "
                + ("(reflecting the confidence sweep). " if cfg.get("sweep") else ". ")
                + "Lower values at higher confidence indicate reads the classifier could not "
                  "place above the threshold, not necessarily absence of signal.")
        blocks.append(_table(["Report", "% classified"], rows, "Percent of reads classified.",
                             "classified"))

    top = _top_taxa(outdir)
    if top:
        blocks.append("\\subsection*{Community composition}")
        rows = [[t.get("name", ""), t.get("sample", ""),
                 t.get("new_est_reads", ""), t.get("fraction_total_reads", "")]
                for t in top]
        blocks.append(_table(["Taxon", "Sample", "Reads", "Fraction"], rows,
                             "Most abundant taxa (Bracken).", "abund", wrap_first=True))
        bp = os.path.join(outdir, "stats", "composition_barplot.png")
        if os.path.isfile(bp):
            blocks.append(_fig(bp, repdir, "Community composition across samples (top taxa).",
                               "barplot"))

    div = _read_json(os.path.join(outdir, "stats", "diversity.json"))
    if div:
        blocks.append("\\subsection*{Diversity}")
        alpha = _alpha_table(outdir)
        if alpha:
            header, ar = alpha
            blocks.append(_table(header, ar,
                                 "Alpha-diversity per sample.", "alpha"))
        expl = div.get("pcoa_explained")
        pc = os.path.join(outdir, "stats", "pcoa.png")
        if expl and os.path.isfile(pc) and os.path.getsize(pc) > 0:
            pcts = ", ".join("{:.1f}\\%".format(x * 100) for x in expl[:2])
            blocks.append(
                "Principal-coordinates analysis of Bray--Curtis dissimilarity explained "
                f"{pcts} of the variance on the first two axes "
                "(Figure~\\ref{fig:pcoa}).")
            blocks.append(_fig(pc, repdir, "PCoA ordination (Bray--Curtis).", "pcoa"))

    da = _read_json(os.path.join(outdir, "stats", "differential_abundance.json"))
    if da:
        sm = da.get("summary", {})
        blocks.append("\\subsection*{Differential abundance}")
        sig = da.get("significant_taxa", [])
        blocks.append(
            f"Comparing {esc(sm.get('group_a','A'))} (n={sm.get('n_a','?')}) with "
            f"{esc(sm.get('group_b','B'))} (n={sm.get('n_b','?')}) over {sm.get('n_taxa','?')} "
            f"taxa, {sm.get('n_significant',0)} taxon(s) were significant at FDR "
            f"{sm.get('fdr','')} ({sm.get('n_permutations','?')}-permutation test on "
            f"centred-log-ratio abundances)."
            + (f" Significant taxa: {esc(', '.join(sig))}." if sig else
               " No taxon passed the FDR threshold, consistent with the available replication."))
        vol = os.path.join(outdir, "stats", "differential_volcano.png")
        if os.path.isfile(vol) and os.path.getsize(vol) > 0:
            blocks.append(_fig(vol, repdir, "Differential-abundance volcano plot.", "volcano"))

    if mods.get("reconcile"):
        for jp in sorted(glob.glob(os.path.join(outdir, "reconcile", "*.reconcile.json"))):
            s = _read_json(jp)
            if not s:
                continue
            c = s.get("taxa_concordance", {})
            blocks.append("\\subsection*{Read--contig reconciliation}")
            blocks.append(
                f"For sample {esc(s.get('sample'))}, {s.get('n_contigs_classified')} of "
                f"{s.get('n_contigs')} contigs were classified; {c.get('both',0)} taxa were "
                f"supported by both reads and contigs, {c.get('reads_only',0)} by reads only, "
                f"and {c.get('contigs_only',0)} by contigs only.")
            break

    if not blocks:
        blocks.append("No result files were found; run the workflow before generating the paper.")
    return "\n\n".join(blocks)


def _discussion(cfg: Dict[str, Any], outdir: str) -> str:
    mods = cfg.get("modules", {})
    p = ["This analysis should be interpreted in light of the methods and their assumptions."]
    if mods.get("classify"):
        p.append(
            "Read-based classification is bounded by the completeness of the reference "
            "database: taxa absent from the database cannot be detected, and database growth "
            "can push assignments toward higher ranks. The reported classified fraction should "
            "therefore be read as a lower bound on community coverage.")
    if cfg.get("sweep"):
        p.append(
            "The confidence sweep makes the precision--recall trade-off explicit: higher "
            "thresholds reduce spurious assignments at the cost of sensitivity. The appropriate "
            "operating point depends on whether the goal is discovery or confident detection.")
    if mods.get("differential"):
        da = _read_json(os.path.join(outdir, "stats", "differential_abundance.json"))
        n = (da or {}).get("summary", {}).get("n_significant", 0)
        if n:
            p.append(
                "The differentially abundant taxa identified here are candidate signals; "
                "because abundances are compositional, the centred-log-ratio framing avoids "
                "spurious negative correlations, but findings warrant validation in an "
                "independent cohort.")
        else:
            p.append(
                "No taxa reached significance. With limited replication a permutation test has "
                "low power (a two-versus-two design cannot reach a small p-value at all), so "
                "this is best read as an absence of evidence rather than evidence of absence; "
                "more biological replicates would be needed to detect moderate effects.")
    if mods.get("assembly"):
        p.append(
            "Assembly-based results depend on sequencing depth and community evenness; "
            "low-abundance or highly similar genomes may remain fragmented, and contig-level "
            "taxonomy reflects identity rather than abundance.")
    p.append(
        "All tool versions, parameters and commands are recorded in the run manifest, so the "
        "analysis is fully reproducible and these limitations can be addressed by targeted "
        "re-runs (deeper sequencing, a more complete database, or additional replicates).")
    return "\n\n".join(p)


# --------------------------------------------------------------------------- #
# Document assembly + compile                                                 #
# --------------------------------------------------------------------------- #
def build_tex(cfg: Dict[str, Any], manifest: Dict[str, Any]) -> str:
    outdir = report._outdir(cfg)
    design = _design(cfg)
    refs = report.citations_for(cfg)
    today = _dt.date.today().isoformat()
    bib = "\n".join(f"\\bibitem{{ref{i}}} {esc(c)}" for i, c in enumerate(refs, 1))

    preamble = (
        "\\documentclass[11pt]{article}\n"
        "\\usepackage[utf8]{inputenc}\n"
        "\\usepackage[margin=1in]{geometry}\n"
        "\\usepackage{graphicx}\n"
        "\\usepackage{hyperref}\n"
        "\\usepackage{float}\n"
        "\\setlength{\\parskip}{0.5em}\n"
        "\\setlength{\\parindent}{0pt}\n"
    )
    head = (
        f"\\title{{{_title(cfg, design)}}}\n"
        f"\\author{{Generated by metagx {esc(__version__)}}}\n"
        f"\\date{{{today}}}\n"
        "\\begin{document}\n\\maketitle\n"
    )
    abstract = ("\\begin{abstract}\n" + _abstract(cfg, design, outdir) + "\n\\end{abstract}\n")
    body = (
        "\\section{Introduction}\n" + _introduction(cfg, design) + "\n\n"
        "\\section{Methods}\n" + _methods(cfg, manifest, design) + "\n\n"
        "\\section{Results}\n" + _results(cfg, manifest, outdir) + "\n\n"
        "\\section{Discussion}\n" + _discussion(cfg, outdir) + "\n\n"
    )
    refs_block = ""
    if bib:
        # \clearpage flushes any still-pending floats before References, so a table/figure can
        # never land after the bibliography (tables/figures also use [H] to pin them in place).
        refs_block = (
            "\\clearpage\n\\begin{thebibliography}{99}\n" + bib + "\n\\end{thebibliography}\n"
        )
    return preamble + head + abstract + body + refs_block + "\\end{document}\n"


def _compile(tex_path: str) -> Optional[str]:
    """Compile a .tex with pdflatex (twice for refs). Returns the pdf path or None."""
    engine = shutil.which("pdflatex")
    if not engine:
        return None
    workdir = os.path.dirname(tex_path)
    name = os.path.basename(tex_path)
    pdf = os.path.splitext(tex_path)[0] + ".pdf"
    for _ in range(2):
        try:
            subprocess.run([engine, "-interaction=nonstopmode", "-halt-on-error", name],
                           cwd=workdir, capture_output=True, text=True, timeout=180)
        except (subprocess.SubprocessError, OSError):
            return None
    # tidy aux files
    base = os.path.splitext(name)[0]
    for ext in (".aux", ".log", ".out", ".toc"):
        f = os.path.join(workdir, base + ext)
        if os.path.exists(f):
            try:
                os.remove(f)
            except OSError:
                pass
    return pdf if os.path.isfile(pdf) else None


def generate(cfg: Dict[str, Any], compile_pdf: bool = True) -> Dict[str, Any]:
    """Write an IMRaD paper (paper.tex) and, if pdflatex is available, paper.pdf.

    Returns ``{"paths": {...}, "pdf": <path or None>, "compiled": bool}``.
    """
    outdir = report._outdir(cfg)
    repdir = os.path.join(outdir, "report")
    os.makedirs(repdir, exist_ok=True)
    manifest = report.build_manifest(cfg)
    tex = build_tex(cfg, manifest)
    tex_path = os.path.join(repdir, "paper.tex")
    with open(tex_path, "w") as fh:
        fh.write(tex)
    paths = {"paper_tex": tex_path}
    pdf = _compile(tex_path) if compile_pdf else None
    if pdf:
        paths["paper_pdf"] = pdf
    return {"paths": paths, "pdf": pdf, "compiled": bool(pdf)}
