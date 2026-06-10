"""Prevalence-based decontamination from negative/blank controls.

A taxon that appears in the negative controls as often as (or more often than) in the real
samples is most likely a reagent/lab contaminant, not biology — the core idea of the decontam
"prevalence" method. We flag such taxa and write a cleaned abundance table with them removed
from the real samples. Pure-Python (csv only); helpers are unit-testable.
"""

import csv
import json
from typing import Dict, List, Optional


def flag_contaminants(rows: List[dict], controls: set,
                      label: Optional[str] = None) -> Dict[str, dict]:
    """Per-taxon prevalence in controls vs real samples + a contaminant flag."""
    if label and any(r.get("label") == label for r in rows):
        rows = [r for r in rows if r.get("label") == label]
    all_samples = {r.get("sample") for r in rows}
    real = all_samples - controls
    n_ctrl, n_real = len(controls), len(real)

    taxa: Dict[str, dict] = {}
    for r in rows:
        frac = float(r.get("fraction_total_reads", 0) or 0)
        if frac <= 0:
            continue
        d = taxa.setdefault(r.get("name"), {"ctrl": set(), "real": set()})
        (d["ctrl"] if r.get("sample") in controls else d["real"]).add(r.get("sample"))

    out: Dict[str, dict] = {}
    for name, d in taxa.items():
        cp = len(d["ctrl"]) / n_ctrl if n_ctrl else 0.0
        sp = len(d["real"]) / n_real if n_real else 0.0
        # contaminant: present in ≥1 control AND not more prevalent among real samples
        out[name] = {"control_prevalence": round(cp, 3), "sample_prevalence": round(sp, 3),
                     "flagged": bool(len(d["ctrl"]) > 0 and cp >= sp)}
    return out


def run(bracken_path: str, controls: List[str], out_flagged: str, out_cleaned: str,
        label: Optional[str] = None) -> dict:
    controls = set(controls)
    with open(bracken_path) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    flags = flag_contaminants(rows, controls, label)
    flagged_taxa = {name for name, d in flags.items() if d["flagged"]}

    with open(out_flagged, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["taxon", "control_prevalence", "sample_prevalence", "flagged"])
        for name in sorted(flags):
            d = flags[name]
            w.writerow([name, d["control_prevalence"], d["sample_prevalence"],
                        "yes" if d["flagged"] else "no"])

    fields = list(rows[0].keys()) if rows else ["sample", "name", "fraction_total_reads"]
    with open(out_cleaned, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for r in rows:
            if r.get("sample") in controls or r.get("name") in flagged_taxa:
                continue
            w.writerow(r)
    return {"n_flagged": len(flagged_taxa), "flagged": sorted(flagged_taxa)}


if "snakemake" in globals():  # pragma: no cover - exercised inside the workflow
    run(
        bracken_path=snakemake.input.bracken,                    # noqa: F821
        controls=list(snakemake.params.controls),                # noqa: F821
        out_flagged=snakemake.output.flagged,                    # noqa: F821
        out_cleaned=snakemake.output.cleaned,                    # noqa: F821
        label=getattr(snakemake.params, "label", None),          # noqa: F821
    )
