"""Side-by-side comparison of the SAME biological sample across sequencing platforms.

Importable engine behind `metagx compare` and `scripts/compare_platforms.py`. Give it
a TSV manifest describing one row per platform (ONT / Illumina / PacBio HiFi / PacBio CLR /
Ion Torrent / ...), each pointing at whatever outputs you have for that platform, and it
produces an integrated comparison: classification accuracy AND assembly contiguity AND
reference recovery, side by side.

Use it whenever you have multi-platform sequencing of one sample (a mock community, a
benchmark, a re-sequenced isolate) and want to see what each technology buys you.

Manifest columns (tab-separated, header required; leave a cell blank to skip that block):

    label            human-readable name, e.g. "ONT (long, noisy)"
    platform         ont | illumina | pacbio_hifi | pacbio_clr   (selects minimap2 preset)
    platform_class   short | long                                (grouping/colour)
    kreport          path to a kraken2 report  -> classification block
    reads            path(s) to input reads (comma-sep for R1,R2) -> read-length + concordance
    contigs          path to assembly FASTA    -> assembly block (may be an empty file)
    reference        path to the ground-truth reference FASTA -> length-bias + breadth

Outputs (results/experiments/cross_platform_comparison/ by default):
    comparison_table.tsv, comparison.json, and fig_*.png

Run with the repo .venv (matplotlib/pandas/numpy) and the bio env on PATH (minimap2/samtools):
    PATH="$HOME/miniconda3/envs/metagx-bio/bin:$PATH" \\
        .venv/bin/python -m metagx.cli compare [--manifest m.tsv] [--outdir d/]
"""
from __future__ import annotations

import csv
import gzip
import json
import math
import os
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass, field
from typing import Dict, List, Optional

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_MANIFEST = os.path.join(ROOT, "config", "cross_platform.manifest.tsv")
DEFAULT_OUTDIR = os.path.join(ROOT, "results", "experiments", "cross_platform_comparison")

PRESET = {"illumina": "sr", "ont": "map-ont", "pacbio_hifi": "map-hifi", "pacbio_clr": "map-pb"}
PALETTE = {"short": "#DD8452", "long": "#4C72B0"}
PER_PLATFORM_COLOR = {
    "ont": "#4C72B0", "illumina": "#DD8452", "pacbio_clr": "#55A868", "pacbio_hifi": "#C44E52",
}


# --------------------------------------------------------------------------- spec
@dataclass
class Spec:
    label: str
    platform: str
    platform_class: str
    kreport: str
    reads: List[str]
    contigs: str
    reference: str


def load_manifest(path: str) -> List[Spec]:
    specs: List[Spec] = []
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if not row.get("label") or row["label"].startswith("#"):
                continue
            reads = [r.strip() for r in (row.get("reads") or "").split(",") if r.strip()]
            specs.append(Spec(
                label=row["label"].strip(),
                platform=(row.get("platform") or "").strip(),
                platform_class=(row.get("platform_class") or "long").strip(),
                kreport=(row.get("kreport") or "").strip(),
                reads=reads,
                contigs=(row.get("contigs") or "").strip(),
                reference=(row.get("reference") or "").strip(),
            ))
    return specs


def _abs(p: str) -> str:
    return p if os.path.isabs(p) else os.path.join(ROOT, p)


# --------------------------------------------------------------------------- parsing
def genome_lengths(ref_path: str) -> Dict[str, int]:
    """Species description (header minus accession) -> bp."""
    lengths: Dict[str, int] = {}
    name: Optional[str] = None
    with open(ref_path) as fh:
        for line in fh:
            if line.startswith(">"):
                desc = line[1:].strip()
                name = desc.split(None, 1)[1] if " " in desc else desc
                lengths[name] = 0
            elif name is not None:
                lengths[name] += len(line.strip())
    return lengths


def genome_lengths_by_acc(ref_path: str) -> Dict[str, int]:
    """Accession token -> bp (minimap2 target names are the accession)."""
    lengths: Dict[str, int] = {}
    acc: Optional[str] = None
    with open(ref_path) as fh:
        for line in fh:
            if line.startswith(">"):
                acc = line[1:].split()[0]
                lengths[acc] = 0
            elif acc is not None:
                lengths[acc] += len(line.strip())
    return lengths


def species_counts(report_path: str) -> pd.Series:
    df = pd.read_csv(report_path, sep="\t", header=None, dtype=str)
    reads = pd.to_numeric(df[1], errors="coerce").fillna(0).astype(int)
    rank = df.iloc[:, -3].str.strip()
    name = df.iloc[:, -1].str.strip()
    out = pd.DataFrame({"name": name, "reads": reads})[rank == "S"]
    return out.groupby("name")["reads"].sum()


def classified_unclassified(report_path: str):
    df = pd.read_csv(report_path, sep="\t", header=None, dtype=str)
    reads = pd.to_numeric(df[1], errors="coerce").fillna(0).astype(int)
    rank = df.iloc[:, -3].str.strip()
    name = df.iloc[:, -1].str.strip().str.lower()
    unclass = int(reads[(rank == "U") | (name == "unclassified")].sum())
    root = int(reads[rank == "R"].sum())
    return root, unclass


def read_lengths(paths: List[str], cap: int = 200_000) -> List[int]:
    lengths: List[int] = []
    for p in paths:
        ap = _abs(p)
        if not os.path.exists(ap):
            continue
        is_fastq = any(ap.endswith(e) for e in (".fastq", ".fq", ".fastq.gz", ".fq.gz"))
        opener = gzip.open if ap.endswith(".gz") else open
        if is_fastq:
            with opener(ap, "rt") as fh:
                for i, line in enumerate(fh):
                    if i % 4 == 1:
                        lengths.append(len(line.strip()))
        else:  # fasta
            cur = 0
            with opener(ap, "rt") as fh:
                for line in fh:
                    if line.startswith(">"):
                        if cur:
                            lengths.append(cur)
                        cur = 0
                    else:
                        cur += len(line.strip())
            if cur:
                lengths.append(cur)
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


def contig_lengths(fa: str) -> List[int]:
    lens, cur = [], 0
    if not fa or not os.path.exists(_abs(fa)):
        return lens
    with open(_abs(fa)) as fh:
        for line in fh:
            if line.startswith(">"):
                if cur:
                    lens.append(cur)
                cur = 0
            else:
                cur += len(line.strip())
    if cur:
        lens.append(cur)
    return lens


def shannon(counts: np.ndarray) -> float:
    p = counts[counts > 0].astype(float)
    if p.sum() == 0:
        return 0.0
    p = p / p.sum()
    return float(-(p * np.log(p)).sum())


def simpson(counts: np.ndarray) -> float:
    p = counts[counts > 0].astype(float)
    if p.sum() == 0:
        return 0.0
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


def _have(tool: str) -> bool:
    return shutil.which(tool) is not None


def reference_breadth(contigs: str, reference: str):
    """minimap2 contigs->reference; return (% of reference bp covered, genomes recovered).

    A genome counts as 'recovered' if >=50% of its length is covered by contigs.
    """
    if not (contigs and os.path.exists(_abs(contigs)) and os.path.getsize(_abs(contigs)) > 0):
        return 0.0, 0
    if not (reference and _have("minimap2")):
        return float("nan"), -1
    glen = genome_lengths_by_acc(_abs(reference))
    covered: Dict[str, List] = {a: [] for a in glen}
    try:
        proc = subprocess.run(
            ["minimap2", "-x", "asm20", "--secondary=no", _abs(reference), _abs(contigs)],
            capture_output=True, text=True, timeout=600,
        )
    except Exception:
        return float("nan"), -1
    for line in proc.stdout.splitlines():
        f = line.split("\t")
        if len(f) < 9:
            continue
        tname, tstart, tend = f[5], int(f[7]), int(f[8])
        if tname in covered:
            covered[tname].append((tstart, tend))
    total_bp = sum(glen.values())
    cov_bp = 0
    recovered = 0
    for acc, ivs in covered.items():
        if not ivs:
            continue
        ivs.sort()
        merged = [list(ivs[0])]
        for s, e in ivs[1:]:
            if s <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], e)
            else:
                merged.append([s, e])
        gcov = sum(e - s for s, e in merged)
        cov_bp += gcov
        if glen[acc] and gcov / glen[acc] >= 0.5:
            recovered += 1
    return round(100 * cov_bp / total_bp, 2) if total_bp else 0.0, recovered


def read_mapping_rate(reads: List[str], contigs: str, platform: str):
    """minimap2 reads->contigs; return % of reads that map (assembly concordance)."""
    if not (contigs and os.path.exists(_abs(contigs)) and os.path.getsize(_abs(contigs)) > 0):
        return 0.0
    if not (reads and _have("minimap2") and _have("samtools")):
        return float("nan")
    preset = PRESET.get(platform, "map-ont")
    rpaths = [_abs(r) for r in reads if os.path.exists(_abs(r))]
    if not rpaths:
        return float("nan")
    try:
        with tempfile.NamedTemporaryFile(suffix=".bam", delete=False) as tmp:
            bam = tmp.name
        mm = subprocess.Popen(
            ["minimap2", "-a", "-x", preset, "-t", "4", _abs(contigs), *rpaths],
            stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
        )
        with open(bam, "wb") as out:
            subprocess.run(["samtools", "view", "-b", "-"], stdin=mm.stdout, stdout=out,
                           stderr=subprocess.DEVNULL, check=True)
        mm.wait()
        fs = subprocess.run(["samtools", "flagstat", bam], capture_output=True, text=True)
        os.unlink(bam)
    except Exception:
        return float("nan")
    total = mapped = 0
    for line in fs.stdout.splitlines():
        if " in total " in line:
            total = int(line.split()[0])
        elif "primary mapped (" in line or (" mapped (" in line and "primary" not in line and total and mapped == 0):
            mapped = int(line.split()[0])
    return round(100 * mapped / total, 2) if total else float("nan")


# --------------------------------------------------------------------------- analyze
@dataclass
class Result:
    spec: Spec
    counts: Optional[pd.Series]
    lengths: List[int]
    clens: List[int]
    row: Dict[str, object] = field(default_factory=dict)


def analyze(spec: Spec) -> Result:
    row: Dict[str, object] = {"platform": spec.label, "platform_class": spec.platform_class}
    counts = None
    lengths = read_lengths(spec.reads) if spec.reads else []

    # classification block
    if spec.kreport and os.path.exists(_abs(spec.kreport)):
        counts = species_counts(_abs(spec.kreport))
        classified, unclassified = classified_unclassified(_abs(spec.kreport))
        total = classified + unclassified
        arr = counts.to_numpy()
        rich = int((arr > 0).sum())
        H = shannon(arr)
        dom = counts.sort_values(ascending=False)
        row.update({
            "reads": total,
            "classified_pct": round(100 * classified / total, 2) if total else 0.0,
            "species_recovered": rich,
            "shannon": round(H, 3),
            "pielou": round(H / math.log(rich), 3) if rich > 1 else 0.0,
            "dominant_taxon": dom.index[0] if len(dom) else "—",
            "dominant_share_pct": round(100 * dom.iloc[0] / classified, 1) if classified and len(dom) else 0.0,
        })
        if spec.reference and os.path.exists(_abs(spec.reference)):
            glen = genome_lengths(_abs(spec.reference))
            gl = np.array([glen[n] for n in glen], dtype=float)
            rd = np.array([int(counts.get(n, 0)) for n in glen], dtype=float)
            rho = spearman(gl, rd)
            row["length_bias_spearman"] = round(rho, 3) if not math.isnan(rho) else None
            row["species_total"] = len(glen)

    # read-length block
    if lengths:
        row.update({
            "read_len_median": int(np.median(lengths)),
            "read_n50": n50(lengths),
            "read_len_max": int(np.max(lengths)),
        })

    # assembly block
    clens = contig_lengths(spec.contigs)
    if spec.contigs:
        if clens:
            breadth, recov = reference_breadth(spec.contigs, spec.reference)
            mrate = read_mapping_rate(spec.reads, spec.contigs, spec.platform)
            row.update({
                "asm_contigs": len(clens),
                "asm_total_bp": sum(clens),
                "asm_n50": n50(clens),
                "asm_longest": max(clens),
                "ref_breadth_pct": breadth,
                "genomes_recovered": recov,
                "read_concordance_pct": mrate,
            })
        else:
            row.update({
                "asm_contigs": 0, "asm_total_bp": 0, "asm_n50": 0, "asm_longest": 0,
                "ref_breadth_pct": 0.0, "genomes_recovered": 0, "read_concordance_pct": 0.0,
                "asm_note": "no assembly (insufficient depth)",
            })
    return Result(spec, counts, lengths, clens, row)


# --------------------------------------------------------------------------- figures
def _short(label: str) -> str:
    return label.split(" (")[0]


def fig_assembly(results: List[Result], path: str) -> None:
    asm = [r for r in results if "asm_contigs" in r.row]
    if not asm:
        return
    labels = [_short(r.spec.label) for r in asm]
    colors = [PER_PLATFORM_COLOR.get(r.spec.platform, "#777") for r in asm]
    n50s = [r.row["asm_n50"] for r in asm]
    ncon = [r.row["asm_contigs"] for r in asm]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.6))
    ax1.bar(labels, n50s, color=colors)
    ax1.set_ylabel("contig N50 (bp)")
    ax1.set_title("Assembly contiguity (higher = better)")
    for i, v in enumerate(n50s):
        ax1.text(i, v, f"{v:,}", ha="center", va="bottom", fontsize=9)
    ax1.tick_params(axis="x", rotation=15)
    ax2.bar(labels, ncon, color=colors)
    ax2.set_ylabel("number of contigs (lower = better)")
    ax2.set_title("Assembly fragmentation")
    for i, v in enumerate(ncon):
        ax2.text(i, v, str(v), ha="center", va="bottom", fontsize=9)
    ax2.tick_params(axis="x", rotation=15)
    fig.suptitle("Long-read vs short-read assembly — same 30-genome reference", y=1.02)
    fig.tight_layout()
    fig.savefig(path, dpi=130, bbox_inches="tight")
    plt.close(fig)


def fig_recovery(results: List[Result], path: str) -> None:
    asm = [r for r in results if "ref_breadth_pct" in r.row]
    if not asm:
        return
    labels = [_short(r.spec.label) for r in asm]
    colors = [PER_PLATFORM_COLOR.get(r.spec.platform, "#777") for r in asm]
    breadth = [r.row["ref_breadth_pct"] for r in asm]
    concord = [r.row.get("read_concordance_pct") or 0 for r in asm]
    x = np.arange(len(labels))
    w = 0.38
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.bar(x - w / 2, breadth, w, label="reference breadth (% of 30 genomes recovered)", color=colors, alpha=0.95)
    ax.bar(x + w / 2, concord, w, label="read concordance (% reads mapping to own contigs)",
           color=colors, alpha=0.5, hatch="//")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=15, ha="right")
    ax.set_ylabel("percent")
    ax.set_ylim(0, 105)
    ax.set_title("Genome reconstruction: reference breadth and read concordance")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def fig_read_length(results: List[Result], path: str) -> None:
    have = [r for r in results if r.lengths]
    if not have:
        return
    fig, ax = plt.subplots(figsize=(9, 5))
    data = [np.array(r.lengths) for r in have]
    labels = [_short(r.spec.label) for r in have]
    parts = ax.violinplot([np.log10(d[d > 0]) for d in data], showmedians=True, widths=0.8)
    for i, pc in enumerate(parts["bodies"]):
        pc.set_facecolor(PER_PLATFORM_COLOR.get(have[i].spec.platform, "#777"))
        pc.set_alpha(0.7)
    ax.set_xticks(range(1, len(labels) + 1))
    ax.set_xticklabels(labels, rotation=15, ha="right")
    yt = [2, 2.5, 3, 3.5, 4, 4.5]
    ax.set_yticks(yt)
    ax.set_yticklabels([f"{int(10**v):,}" for v in yt])
    ax.set_ylabel("read length (bp, log)")
    ax.set_title("Read-length distribution by platform")
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


# --------------------------------------------------------------------------- main
def run(manifest: Optional[str] = None, outdir: Optional[str] = None) -> str:
    """Run the comparison; return the output directory. Engine used by CLI + script."""
    manifest = manifest or DEFAULT_MANIFEST
    outdir = outdir or DEFAULT_OUTDIR
    os.makedirs(outdir, exist_ok=True)

    specs = load_manifest(manifest)
    results = [analyze(s) for s in specs]
    table = pd.DataFrame([r.row for r in results])

    tsv = os.path.join(outdir, "comparison_table.tsv")
    table.to_csv(tsv, sep="\t", index=False)
    with open(os.path.join(outdir, "comparison.json"), "w") as fh:
        json.dump({"manifest": manifest, "platforms": [r.row for r in results]}, fh, indent=2)

    fig_read_length(results, os.path.join(outdir, "fig_read_length.png"))
    fig_assembly(results, os.path.join(outdir, "fig_assembly.png"))
    fig_recovery(results, os.path.join(outdir, "fig_recovery.png"))

    with pd.option_context("display.max_columns", None, "display.width", 200):
        print(table.to_string(index=False))
    print(f"\nWrote {tsv} + comparison.json + fig_*.png to {outdir}")
    return outdir


def main(argv: Optional[List[str]] = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    manifest = argv[0] if len(argv) > 0 else None
    outdir = argv[1] if len(argv) > 1 else None
    run(manifest, outdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
