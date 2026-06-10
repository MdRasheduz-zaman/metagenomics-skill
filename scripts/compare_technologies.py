#!/usr/bin/env python3
"""Cross-sequencing-technology comparison on a shared 30-genome viral reference.

NOTE: `scripts/compare_platforms.py` is the generalized, manifest-driven successor to this
script — it covers classification AND assembly for an arbitrary set of platforms via a TSV
manifest. This script remains as the focused classification-only experiment 07.


All four datasets (ONT, Illumina, PacBio CLR, PacBio HiFi) were simulated from the
*same* `data/genomes.fasta`, so differences in the kraken2 profiles are attributable to
the sequencing technology / simulator, not the underlying community. This script
quantifies those differences and writes a table + figures + an IMRaD markdown report to
`results/experiments/cross_technology_comparison/`.

Metrics per technology:
  - input read count, read-length distribution (min/median/mean/N50/max)
  - classified vs unclassified fraction (kraken2 confidence 0.0)
  - species recovered out of the 30 reference genomes (sensitivity)
  - Shannon / Simpson / Pielou evenness of the species profile
  - dominant taxon and its share
  - genome-length bias: Spearman correlation between per-genome assigned reads and
    reference genome length (uniform-coverage simulators inflate the longest genome)

Run with the repo's .venv python (matplotlib/pandas/numpy live there):
    .venv/bin/python scripts/compare_technologies.py
"""
from __future__ import annotations

import gzip
import json
import math
import os
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REF = os.path.join(ROOT, "data", "genomes.fasta")
OUTDIR = os.path.join(ROOT, "results", "experiments", "cross_technology_comparison")


@dataclass
class Tech:
    key: str
    label: str
    kreport: str
    reads: List[str]  # input read files (fasta or fastq.gz); R1 only for paired is fine for length
    read_format: str  # "fasta" | "fastq"
    simulator: str
    platform_class: str  # "short" | "long"


TECHS: List[Tech] = [
    Tech(
        key="ont_sim",
        label="ONT (long, noisy)",
        kreport="results/ont_sim/kraken2/ont_sim.confidence_0.0.kreport",
        reads=["data/simulated_metagenomic_reads.fasta"],
        read_format="fasta",
        simulator="bundled ONT sim (per-genome read count)",
        platform_class="long",
    ),
    Tech(
        key="illumina_sim",
        label="Illumina (short, accurate)",
        kreport="results/illumina_sim/kraken2/illumina_sim.confidence_0.0.kreport",
        reads=["data/illumina_sim/illumina_sim_R1.fastq.gz"],
        read_format="fastq",
        simulator="wgsim (uniform coverage)",
        platform_class="short",
    ),
    Tech(
        key="pacbio_clr_sim",
        label="PacBio CLR (long, noisy)",
        kreport="results/pacbio_clr_sim/kraken2/pacbio_clr_sim.confidence_0.0.kreport",
        reads=["data/pacbio_sim/pacbio_clr_simlord.fastq.gz"],
        read_format="fastq",
        simulator="SimLoRD 1-pass (uniform coverage)",
        platform_class="long",
    ),
    Tech(
        key="pacbio_hifi_sim",
        label="PacBio HiFi (long, accurate)",
        kreport="results/pacbio_hifi_sim/kraken2/pacbio_hifi_sim.confidence_0.0.kreport",
        reads=["data/pacbio_sim/pacbio_hifi.fastq.gz"],
        read_format="fastq",
        simulator="SimLoRD multipass (uniform coverage)",
        platform_class="long",
    ),
]


# --------------------------------------------------------------------------- parsing
def genome_lengths(ref_path: str) -> Dict[str, int]:
    """Map species description (header minus the accession token) -> total bp.

    kraken2 species names match the FASTA description after the accession, e.g.
    '>NC_010356.1 Glossina pallidipes ... genome' -> 'Glossina pallidipes ... genome'.
    """
    lengths: Dict[str, int] = {}
    name: Optional[str] = None
    n = 0
    with open(ref_path) as fh:
        for line in fh:
            if line.startswith(">"):
                # drop '>' and the first whitespace-delimited token (accession)
                desc = line[1:].strip()
                name = desc.split(None, 1)[1] if " " in desc else desc
                lengths[name] = 0
            elif name is not None:
                lengths[name] += len(line.strip())
                n += 1
    return lengths


def species_counts(report_path: str) -> pd.Series:
    """Species-level reads-in-clade, robust to --report-minimizer-data extra columns."""
    df = pd.read_csv(report_path, sep="\t", header=None, dtype=str)
    reads = pd.to_numeric(df[1], errors="coerce").fillna(0).astype(int)
    rank = df.iloc[:, -3].str.strip()
    name = df.iloc[:, -1].str.strip()
    out = pd.DataFrame({"name": name, "reads": reads})[rank == "S"]
    return out.groupby("name")["reads"].sum()


def classified_unclassified(report_path: str) -> Tuple[int, int]:
    """Return (classified_reads, unclassified_reads) from the U and R rows."""
    df = pd.read_csv(report_path, sep="\t", header=None, dtype=str)
    reads = pd.to_numeric(df[1], errors="coerce").fillna(0).astype(int)
    rank = df.iloc[:, -3].str.strip()
    name = df.iloc[:, -1].str.strip().str.lower()
    unclass = int(reads[(rank == "U") | (name == "unclassified")].sum())
    root = int(reads[rank == "R"].sum())
    return root, unclass


def read_lengths(paths: List[str], fmt: str, cap: int = 200_000) -> List[int]:
    lengths: List[int] = []
    for p in paths:
        ap = os.path.join(ROOT, p) if not os.path.isabs(p) else p
        if not os.path.exists(ap):
            continue
        if fmt == "fasta":
            cur = 0
            with open(ap) as fh:
                for line in fh:
                    if line.startswith(">"):
                        if cur:
                            lengths.append(cur)
                        cur = 0
                    else:
                        cur += len(line.strip())
            if cur:
                lengths.append(cur)
        else:  # fastq.gz
            opener = gzip.open if ap.endswith(".gz") else open
            with opener(ap, "rt") as fh:
                for i, line in enumerate(fh):
                    if i % 4 == 1:
                        lengths.append(len(line.strip()))
        if len(lengths) > cap:
            break
    return lengths


# --------------------------------------------------------------------------- metrics
def n50(values: List[int]) -> int:
    if not values:
        return 0
    s = sorted(values, reverse=True)
    half = sum(s) / 2
    run = 0
    for v in s:
        run += v
        if run >= half:
            return v
    return s[-1]


def shannon(counts: np.ndarray) -> float:
    p = counts[counts > 0].astype(float)
    p = p / p.sum()
    return float(-(p * np.log(p)).sum())


def simpson(counts: np.ndarray) -> float:
    p = counts[counts > 0].astype(float)
    p = p / p.sum()
    return float(1.0 - (p * p).sum())


def spearman(x: np.ndarray, y: np.ndarray) -> float:
    if len(x) < 3:
        return float("nan")
    rx = pd.Series(x).rank().to_numpy()
    ry = pd.Series(y).rank().to_numpy()
    if rx.std() == 0 or ry.std() == 0:
        return float("nan")
    return float(np.corrcoef(rx, ry)[0, 1])


@dataclass
class Result:
    tech: Tech
    counts: pd.Series
    lengths: List[int]
    classified: int
    unclassified: int
    glengths: Dict[str, int]
    row: Dict[str, object] = field(default_factory=dict)


def analyze(tech: Tech, glengths: Dict[str, int]) -> Result:
    kpath = os.path.join(ROOT, tech.kreport)
    counts = species_counts(kpath)
    classified, unclassified = classified_unclassified(kpath)
    lengths = read_lengths(tech.reads, tech.read_format)
    total = classified + unclassified
    n_species = int((counts > 0).sum())
    arr = counts.to_numpy()
    H = shannon(arr)
    rich = int((arr > 0).sum())
    pielou = H / math.log(rich) if rich > 1 else 0.0
    dom = counts.sort_values(ascending=False)
    dom_name = dom.index[0] if len(dom) else "—"
    dom_share = (dom.iloc[0] / classified * 100) if classified else 0.0

    # genome-length bias: per-genome assigned reads vs reference length
    paired = [(glengths[name], int(counts.get(name, 0))) for name in glengths]
    gl = np.array([a for a, _ in paired], dtype=float)
    rd = np.array([b for _, b in paired], dtype=float)
    rho = spearman(gl, rd)

    row = {
        "technology": tech.label,
        "simulator": tech.simulator,
        "platform_class": tech.platform_class,
        "total_reads": total,
        "classified_reads": classified,
        "classified_pct": round(100 * classified / total, 2) if total else 0.0,
        "species_recovered": n_species,
        "species_total": len(glengths),
        "sensitivity_pct": round(100 * n_species / len(glengths), 1),
        "read_len_median": int(np.median(lengths)) if lengths else 0,
        "read_len_mean": int(np.mean(lengths)) if lengths else 0,
        "read_n50": n50(lengths),
        "read_len_max": int(np.max(lengths)) if lengths else 0,
        "shannon": round(H, 3),
        "simpson": round(simpson(arr), 3),
        "pielou_evenness": round(pielou, 3),
        "dominant_taxon": dom_name,
        "dominant_share_pct": round(dom_share, 1),
        "length_bias_spearman": round(rho, 3) if not math.isnan(rho) else None,
    }
    return Result(tech, counts, lengths, classified, unclassified, glengths, row)


# --------------------------------------------------------------------------- figures
PALETTE = {
    "ont_sim": "#4C72B0",
    "illumina_sim": "#DD8452",
    "pacbio_clr_sim": "#55A868",
    "pacbio_hifi_sim": "#C44E52",
}


def fig_read_length(results: List[Result], path: str) -> None:
    fig, ax = plt.subplots(figsize=(9, 5))
    data = [np.array(r.lengths) for r in results]
    labels = [r.tech.label.split(" (")[0] for r in results]
    parts = ax.violinplot(
        [np.log10(d[d > 0]) for d in data], showmedians=True, widths=0.8
    )
    for i, pc in enumerate(parts["bodies"]):
        pc.set_facecolor(PALETTE[results[i].tech.key])
        pc.set_alpha(0.7)
    ax.set_xticks(range(1, len(labels) + 1))
    ax.set_xticklabels(labels, rotation=15, ha="right")
    ax.set_ylabel("read length (log10 bp)")
    ax.set_title("Read-length distribution by sequencing technology\n(same 30-genome viral reference)")
    yt = [2, 2.5, 3, 3.5, 4, 4.5]
    ax.set_yticks(yt)
    ax.set_yticklabels([f"{int(10**v):,}" for v in yt])
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def fig_classified_and_sensitivity(results: List[Result], path: str) -> None:
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.6))
    labels = [r.tech.label.split(" (")[0] for r in results]
    colors = [PALETTE[r.tech.key] for r in results]
    cl = [r.row["classified_pct"] for r in results]
    sens = [r.row["species_recovered"] for r in results]
    total = results[0].row["species_total"]

    ax1.bar(labels, cl, color=colors)
    ax1.set_ylim(0, 105)
    ax1.set_ylabel("reads classified (%)")
    ax1.set_title("Classified fraction (kraken2, confidence 0.0)")
    for i, v in enumerate(cl):
        ax1.text(i, v + 1.5, f"{v:.1f}%", ha="center", fontsize=9)
    ax1.tick_params(axis="x", rotation=15)

    ax2.bar(labels, sens, color=colors)
    ax2.axhline(total, ls="--", color="grey", lw=1)
    ax2.set_ylim(0, total + 3)
    ax2.set_ylabel(f"species recovered (of {total})")
    ax2.set_title("Detection sensitivity")
    for i, v in enumerate(sens):
        ax2.text(i, v + 0.3, str(v), ha="center", fontsize=9)
    ax2.tick_params(axis="x", rotation=15)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def fig_length_bias(results: List[Result], glengths: Dict[str, int], path: str) -> None:
    fig, axes = plt.subplots(1, len(results), figsize=(4.4 * len(results), 4.4), sharey=False)
    if len(results) == 1:
        axes = [axes]
    names = list(glengths.keys())
    gl = np.array([glengths[n] for n in names], dtype=float)
    for ax, r in zip(axes, results):
        rd = np.array([int(r.counts.get(n, 0)) for n in names], dtype=float)
        ax.scatter(gl, rd + 0.5, s=22, color=PALETTE[r.tech.key], alpha=0.8, edgecolor="k", linewidth=0.3)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("genome length (bp)")
        rho = r.row["length_bias_spearman"]
        ax.set_title(f"{r.tech.label.split(' (')[0]}\nSpearman ρ = {rho}")
        ax.grid(alpha=0.3)
    axes[0].set_ylabel("reads assigned (+0.5, log)")
    fig.suptitle("Genome-length bias: assigned reads vs reference length", y=1.02)
    fig.tight_layout()
    fig.savefig(path, dpi=130, bbox_inches="tight")
    plt.close(fig)


def fig_abundance_heatmap(results: List[Result], path: str) -> None:
    # relative abundance (% of classified) per species x technology
    mat = pd.DataFrame({r.tech.label.split(" (")[0]: r.counts for r in results}).fillna(0.0)
    rel = mat.div(mat.sum(axis=0), axis=1) * 100.0
    rel = rel.loc[rel.max(axis=1).sort_values(ascending=False).index]
    # shorten species labels
    short = [n.replace(", complete genome", "").replace(", complete sequence", "")[:42] for n in rel.index]
    fig, ax = plt.subplots(figsize=(7.5, max(5, 0.32 * len(rel))))
    im = ax.imshow(rel.to_numpy(), aspect="auto", cmap="viridis")
    ax.set_xticks(range(rel.shape[1]))
    ax.set_xticklabels(rel.columns, rotation=20, ha="right")
    ax.set_yticks(range(len(short)))
    ax.set_yticklabels(short, fontsize=7)
    ax.set_title("Relative abundance (% of classified reads)")
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("% of classified")
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


# --------------------------------------------------------------------------- main
def main() -> None:
    os.makedirs(OUTDIR, exist_ok=True)
    glengths = genome_lengths(REF)
    results = [analyze(t, glengths) for t in TECHS]

    table = pd.DataFrame([r.row for r in results])
    tsv = os.path.join(OUTDIR, "comparison_table.tsv")
    table.to_csv(tsv, sep="\t", index=False)

    payload = {
        "reference": "data/genomes.fasta",
        "n_reference_genomes": len(glengths),
        "longest_genome": max(glengths, key=glengths.get),
        "longest_genome_bp": max(glengths.values()),
        "total_reference_bp": sum(glengths.values()),
        "kraken2_confidence": 0.0,
        "technologies": [r.row for r in results],
    }
    with open(os.path.join(OUTDIR, "comparison.json"), "w") as fh:
        json.dump(payload, fh, indent=2)

    fig_read_length(results, os.path.join(OUTDIR, "fig_read_length.png"))
    fig_classified_and_sensitivity(results, os.path.join(OUTDIR, "fig_classified_sensitivity.png"))
    fig_length_bias(results, glengths, os.path.join(OUTDIR, "fig_length_bias.png"))
    fig_abundance_heatmap(results, os.path.join(OUTDIR, "fig_abundance_heatmap.png"))

    print(table.to_string(index=False))
    print(f"\nWrote: {tsv}")
    print(f"Wrote: {os.path.join(OUTDIR, 'comparison.json')}")
    print("Figures: fig_read_length, fig_classified_sensitivity, fig_length_bias, fig_abundance_heatmap")


if __name__ == "__main__":
    main()
