"""Cross-check kraken2 (nucleotide k-mer) against a second, independent classifier
(MetaPhlAn markers or Kaiju protein) at the species level.

Agreement between two orthogonal methods is a strong confidence signal; taxa seen by only
one method flag likely database-completeness false positives (kraken2-only) or divergent
organisms kraken2 missed (consensus-only). Emits a per-sample concordance JSON.

Runs as a Snakemake `script:` (reads snakemake.input/output/params); the parsing helpers
are pure functions so they can be unit-tested directly.
"""

from __future__ import annotations

import json
import re
from typing import Dict, Set


def _norm(name: str) -> str:
    """Normalize a taxon label to a comparable 'genus species' key."""
    s = re.sub(r"^[a-z]__", "", name.strip())      # drop metaphlan s__/g__ prefixes
    s = s.replace("_", " ").strip().lower()
    parts = [p for p in s.split() if p]
    return " ".join(parts[:2])                       # genus + species


def parse_kraken_species(path: str) -> Dict[str, float]:
    """kraken2 kreport species rows (rank 'S') -> {species: percent}."""
    out: Dict[str, float] = {}
    with open(path) as fh:
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if len(cols) >= 6 and cols[3].strip() == "S":
                key = _norm(cols[5])
                if key:
                    out[key] = out.get(key, 0.0) + float(cols[0])
    return out


def parse_metaphlan(path: str) -> Dict[str, float]:
    """MetaPhlAn profile -> {species: relative_abundance%} for s__ rows."""
    out: Dict[str, float] = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            clade = cols[0]
            if "|s__" in clade and "|t__" not in clade:
                sp = clade.split("|s__")[-1]
                try:
                    out[_norm(sp)] = float(cols[2])
                except (IndexError, ValueError):
                    continue
    return out


def parse_kaiju_table(path: str) -> Dict[str, float]:
    """kaiju2table output (file, percent, reads, taxon_id, taxon_name) -> {species: percent}."""
    out: Dict[str, float] = {}
    with open(path) as fh:
        header = fh.readline().lower()
        cols_h = header.rstrip("\n").split("\t")
        try:
            pi = cols_h.index("percent")
        except ValueError:
            pi = 1
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if len(cols) <= pi:
                continue
            name = cols[-1].strip()
            low = name.lower()
            # Drop kaiju's non-species rows. The "cannot be assigned" text is DB-dependent
            # (e.g. "...to a species" vs "...to a (non-viral) species"), so prefix-match it.
            if not name or low == "unclassified" or low.startswith("cannot be assigned"):
                continue
            try:
                val = float(cols[pi])
            except ValueError:
                continue
            # zero-abundance rows are higher-rank catch-alls (e.g. taxid 10239 "Viruses" at
            # species rank with 0 reads), not detected species.
            if val <= 0:
                continue
            out[_norm(name)] = out.get(_norm(name), 0.0) + val
    return out


def concordance(kraken: Dict[str, float], other: Dict[str, float],
                other_name: str, sample: str, top_n: int = 10) -> Dict:
    ks: Set[str] = {k for k in kraken if k}
    os_: Set[str] = {k for k in other if k}
    shared = ks & os_
    union = ks | os_
    top_k = {k for k, _ in sorted(kraken.items(), key=lambda kv: -kv[1])[:top_n]}
    top_o = {k for k, _ in sorted(other.items(), key=lambda kv: -kv[1])[:top_n]}
    return {
        "sample": sample,
        "second_classifier": other_name,
        "n_species_kraken2": len(ks),
        f"n_species_{other_name}": len(os_),
        "n_shared": len(shared),
        "jaccard": round(len(shared) / len(union), 4) if union else 0.0,
        "top_overlap": sorted(top_k & top_o),
        "kraken2_only": sorted(ks - os_)[:50],
        f"{other_name}_only": sorted(os_ - ks)[:50],
    }


def run(kraken_report: str, second_profile: str, classifier: str,
        sample: str, out_json: str) -> Dict:
    kraken = parse_kraken_species(kraken_report)
    if classifier == "metaphlan":
        other = parse_metaphlan(second_profile)
    else:
        other = parse_kaiju_table(second_profile)
    result = concordance(kraken, other, classifier, sample)
    with open(out_json, "w") as fh:
        json.dump(result, fh, indent=2)
    return result


if "snakemake" in globals():  # pragma: no cover - exercised inside the workflow
    run(
        kraken_report=snakemake.input.kraken,            # noqa: F821
        second_profile=snakemake.input.second,           # noqa: F821
        classifier=snakemake.params.classifier,          # noqa: F821
        sample=snakemake.params.sample,                  # noqa: F821
        out_json=snakemake.output.json,                  # noqa: F821
    )
