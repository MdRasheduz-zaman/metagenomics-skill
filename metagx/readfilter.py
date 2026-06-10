"""Taxonomic read filtering before assembly — no external tool (uses our kraken2 per-read
output directly).

Two modes:
  exclude  (default) — drop reads assigned to the given taxids (e.g. host); KEEP everything
                       else including unclassified. This is depletion: it removes a dominant
                       /contaminant fraction so the rest assembles better.
  include            — keep only reads assigned to the given taxids (targeted recovery).
                       Set keep_unclassified=True to also keep unclassified reads — important,
                       since novel/divergent target reads are often unclassified and dropping
                       them biases discovery and fragments the target assembly.

include_children expands each taxid to its whole clade via nodes.dmp (e.g. Viruses 10239).
Abundance (Bracken) should always be computed on the UNfiltered reads; filtering is for
assembly only.
"""

from __future__ import annotations

import gzip
from typing import Dict, List

from .formats import is_gzipped, read_format
from .subsample import _iter_fasta, _iter_fastq


def _open(path: str):
    return gzip.open(path, "rt") if is_gzipped(path) else open(path, "rt")


def parse_read_taxa(kraken_out: str) -> Dict[str, int]:
    d = {}
    with open(kraken_out) as fh:
        for line in fh:
            c = line.split("\t")
            if len(c) < 3:
                continue
            try:
                d[c[1]] = int(c[2])
            except ValueError:
                d[c[1]] = 0
    return d


def _children_map(nodes_path: str) -> Dict[int, List[int]]:
    ch: Dict[int, List[int]] = {}
    try:
        with open(nodes_path) as fh:
            for line in fh:
                p = [x.strip() for x in line.split("|")]
                if len(p) >= 2 and p[0].isdigit():
                    child, par = int(p[0]), int(p[1])
                    if child != par:
                        ch.setdefault(par, []).append(child)
    except OSError:
        pass
    return ch


def expand_taxids(taxids, nodes_path: str) -> set:
    ch = _children_map(nodes_path)
    out, stack = set(taxids), list(taxids)
    while stack:
        for c in ch.get(stack.pop(), []):
            if c not in out:
                out.add(c)
                stack.append(c)
    return out


def _keep(taxid: int, target: set, mode: str, keep_unclassified: bool) -> bool:
    in_target = taxid in target
    if mode == "include":
        return in_target or (keep_unclassified and taxid == 0)
    return not in_target  # exclude: unclassified (0) kept unless 0 is explicitly targeted


def filter_reads(reads: str, kraken_out: str, out: str, taxids, nodes_path: str = "",
                 include_children: bool = True, mode: str = "exclude",
                 keep_unclassified: bool = True) -> dict:
    """Write a filtered copy of ``reads`` (FASTA/FASTQ, ±gz). Returns {total, kept, mode}."""
    read_tax = parse_read_taxa(kraken_out)
    target = set(taxids or [])
    if include_children and nodes_path and target:
        target = expand_taxids(target, nodes_path)
    fmt = read_format(reads)
    total = kept = 0
    with _open(reads) as fin, open(out, "w") as fout:
        if fmt == "fastq":
            for h, s, plus, q in _iter_fastq(fin):
                total += 1
                if _keep(read_tax.get(h[1:].split()[0], 0), target, mode, keep_unclassified):
                    fout.write(f"{h}\n{s}\n{plus}\n{q}\n")
                    kept += 1
        else:
            for header, seq in _iter_fasta(fin):
                total += 1
                if _keep(read_tax.get(header[1:].split()[0], 0), target, mode, keep_unclassified):
                    fout.write(f"{header}\n{seq}\n")
                    kept += 1
    return {"total": total, "kept": kept, "removed": total - kept,
            "mode": mode, "n_target_taxa": len(target)}
