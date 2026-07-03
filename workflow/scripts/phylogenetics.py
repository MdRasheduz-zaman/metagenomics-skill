"""Phylogenetics: MAFFT alignment → optional TrimAl → IQ-TREE 2 or FastTree → stats + plot.

Follows the standard pipeline from the phylogenetics skill (MAFFT / IQ-TREE / FastTree).
ETE3 is optional; basic Newick stats and a simple tree figure use matplotlib when ETE3
is not installed.
"""

import json
import os
import re
import shutil
import subprocess
import sys

# NOTE: no `from __future__ import annotations` — Snakemake's `script:` directive prepends a
# preamble to this file at runtime, which makes a __future__ import no longer the first
# statement (SyntaxError). Plain annotations are fine on Python >=3.10.

# Snakemake script: repo root on path for metagx.registries
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

from metagx import registry  # noqa: E402


def _iqtree_binary() -> str:
    """Resolve the IQ-TREE executable across versions.

    IQ-TREE 2 ships ``iqtree2``; IQ-TREE 3 ships ``iqtree3`` (and sometimes ``iqtree``);
    older/distro builds use ``iqtree``. Hardcoding ``iqtree2`` silently breaks the default
    phylogenetics path the moment the installed IQ-TREE is a different major version — exactly
    the kind of tool-version drift this pipeline must not fall over on. Prefer v2 (the
    registry's reference), then v3, then the generic name.
    """
    for cand in ("iqtree2", "iqtree3", "iqtree"):
        if shutil.which(cand):
            return cand
    raise FileNotFoundError(
        "no IQ-TREE executable found (looked for iqtree2 / iqtree3 / iqtree). "
        "Install IQ-TREE or set phylogenetics.method=fasttree."
    )

MAFFT_METHODS = {
    "auto": ["--auto"],
    "linsi": ["--localpair", "--maxiterate", "1000"],
    "einsi": ["--genafpair", "--maxiterate", "1000"],
    "fftnsi": ["--fftnsi"],
    "fftns": ["--fftns"],
    "retree2": ["--retree", "2"],
}


def count_sequences(fasta_path: str) -> int:
    n = 0
    with open(fasta_path) as fh:
        for line in fh:
            if line.startswith(">"):
                n += 1
    return n


def run_mafft(input_fasta: str, output_fasta: str, mafft_cfg: dict, threads: int) -> None:
    method = str(mafft_cfg.get("method", "auto")).lower()
    extra = list(MAFFT_METHODS.get(method, MAFFT_METHODS["auto"]))
    maxit = int(mafft_cfg.get("maxiterate") or 0)
    if maxit > 0 and method in ("linsi", "einsi"):
        extra = extra[:-1] + [str(maxit)] if extra[-1].isdigit() else extra + ["--maxiterate", str(maxit)]
    cmd = ["mafft", "--thread", str(threads), "--inputorder"] + extra + [input_fasta]
    with open(output_fasta, "w") as out:
        subprocess.run(cmd, check=True, stdout=out, stderr=subprocess.PIPE, text=True)


def run_trimal(input_fasta: str, output_fasta: str, method: str = "automated1") -> bool:
    cmd = ["trimal", f"-{method}", "-in", input_fasta, "-out", output_fasta, "-fasta"]
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        return True
    except (subprocess.CalledProcessError, OSError):
        # trimming is OPTIONAL — fall back to the untrimmed alignment on ANY failure to run trimal
        # (absent => FileNotFoundError, present-but-not-executable => PermissionError, etc.). OSError
        # is the shared parent, so a broken/missing trimal degrades gracefully instead of crashing.
        shutil.copy(input_fasta, output_fasta)
        return False


def run_iqtree(aligned: str, prefix: str, iqtree_cfg: dict, threads: int) -> str:
    managed = {
        "alignment": aligned,
        "prefix": prefix,
        "threads": threads,
        "redo": True,
    }
    args = registry.render_args("iqtree", iqtree_cfg or {}, managed=managed)
    subprocess.run([_iqtree_binary()] + args, check=True)
    tree = f"{prefix}.treefile"
    if not os.path.isfile(tree):
        raise RuntimeError("IQ-TREE did not produce .treefile")
    return tree


def run_fasttree(aligned: str, output_tree: str, ft_cfg: dict) -> None:
    seq_type = str(ft_cfg.get("sequence_type", "nt")).lower()
    model = str(ft_cfg.get("model", "gtr")).lower()
    cmd = ["FastTree"]
    if seq_type == "nt":
        cmd += ["-nt", "-gtr"] if model == "gtr" else ["-nt"]
    else:
        cmd += [f"-{model}" if model in ("lg", "wag", "jtt") else "-lg"]
    cmd.append(aligned)
    with open(output_tree, "w") as out:
        subprocess.run(cmd, check=True, stdout=out, stderr=subprocess.PIPE, text=True)


def newick_stats(tree_path: str) -> dict:
    with open(tree_path) as fh:
        text = fh.read().strip()
    # Leaves are names preceded by '(' or ',' — NOT internal-node support values,
    # which follow ')'. (FastTree writes support like ")0.993:", which a naive
    # "token before a colon" count would wrongly treat as tips.)
    leaves = len(re.findall(r"[(,]\s*([^(),:;\s]+)\s*:", text)) or text.count(",") + 1
    branch_lengths = [float(x) for x in re.findall(r":([0-9.eE+-]+)", text)]
    return {
        "n_leaves": leaves,
        "total_branch_length": round(sum(branch_lengths), 6) if branch_lengths else 0.0,
        "max_branch_length": round(max(branch_lengths), 6) if branch_lengths else 0.0,
    }


def plot_tree(tree_path: str, output_png: str) -> bool:
    try:
        from ete3 import Tree, TreeStyle  # noqa: WPS433

        t = Tree(tree_path)
        t.set_outgroup(t.get_midpoint_outgroup())
        ts = TreeStyle()
        ts.show_leaf_name = True
        ts.mode = "r"
        t.render(output_png, tree_style=ts, w=800, units="px")
        return True
    except Exception:
        pass
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt  # noqa: WPS433

        with open(tree_path) as fh:
            n = fh.read().count(",") + 1
        fig, ax = plt.subplots(figsize=(8, max(3, n * 0.15)))
        ax.text(0.5, 0.5, f"Phylogenetic tree ({n} tips)\nSee {os.path.basename(tree_path)}",
                ha="center", va="center", fontsize=11)
        ax.axis("off")
        fig.savefig(output_png, dpi=120, bbox_inches="tight")
        plt.close(fig)
        return True
    except Exception:
        open(output_png, "a").close()
        return False


def main(
    input_fasta,
    aligned_out,
    tree_out,
    json_out,
    plot_out,
    phylo_cfg,
    mafft_cfg,
    iqtree_cfg,
    fasttree_cfg,
    threads,
):
    os.makedirs(os.path.dirname(aligned_out), exist_ok=True)
    method = str(phylo_cfg.get("method", "iqtree")).lower()
    seq_type = str(phylo_cfg.get("sequence_type", "nt")).lower()
    trim = bool(phylo_cfg.get("trim", True))
    aligned_input = phylo_cfg.get("aligned_input")

    if aligned_input and os.path.isfile(aligned_input):
        aligned = aligned_input
        shutil.copy(aligned, aligned_out)
        aligned_step = "skipped (pre-aligned)"
    else:
        run_mafft(input_fasta, aligned_out, mafft_cfg or {}, threads)
        aligned = aligned_out
        aligned_step = "mafft"

    n_seq = count_sequences(aligned)
    if trim:
        trimmed = aligned_out.replace(".fasta", ".trimmed.fasta")
        trimal_ok = run_trimal(aligned, trimmed, method=str(phylo_cfg.get("trimal_method", "automated1")))
        if trimal_ok:
            aligned = trimmed
            shutil.copy(trimmed, aligned_out)

    use_fasttree = method == "fasttree" or (
        method == "auto" and n_seq > int(phylo_cfg.get("fasttree_threshold", 500))
    )
    prefix = os.path.splitext(tree_out)[0]

    if use_fasttree:
        ft_cfg = dict(fasttree_cfg or {})
        ft_cfg.setdefault("sequence_type", seq_type)
        run_fasttree(aligned, tree_out, ft_cfg)
        tree_method = "fasttree"
    else:
        iq_cfg = dict(iqtree_cfg or {})
        if seq_type == "aa" and iq_cfg.get("model") in (None, "TEST"):
            iq_cfg["model"] = "TEST"
        run_iqtree(aligned, prefix, iq_cfg, threads)
        shutil.copy(f"{prefix}.treefile", tree_out)
        tree_method = _iqtree_binary()   # record the actual IQ-TREE used (iqtree2/3/iqtree)

    stats = newick_stats(tree_out)
    plotted = plot_tree(tree_out, plot_out)

    payload = {
        "n_sequences": n_seq,
        "alignment_step": aligned_step,
        "trimmed": trim,
        "tree_method": tree_method,
        "sequence_type": seq_type,
        "tree_stats": stats,
        "plot_generated": plotted,
        "outputs": {
            "aligned": aligned_out,
            "tree": tree_out,
            "plot": plot_out,
        },
    }
    with open(json_out, "w") as fh:
        json.dump(payload, fh, indent=2)


if __name__ == "__main__":
    sm = snakemake  # noqa: F821
    cfg = sm.params.phylo or {}
    main(
        sm.input.fasta,
        sm.output.aligned,
        sm.output.tree,
        sm.output.json,
        sm.output.plot,
        cfg,
        sm.params.mafft,
        sm.params.iqtree,
        sm.params.fasttree,
        int(sm.params.threads),
    )
