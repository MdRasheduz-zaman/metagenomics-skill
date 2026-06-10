"""Authenticate ancient DNA from mapDamage2 deamination frequency tables.

Authentic ancient DNA carries a post-mortem damage signature: elevated C→T substitutions at
5' read ends and G→A at 3' ends (cytosine deamination). This reads mapDamage's terminal-
position frequencies and emits a verdict; modern contamination shows little/no terminal
excess. Runs as a Snakemake `script:`; the helpers are pure functions for unit testing.
"""

import json


def first_position_freq(path: str) -> float:
    """Return the substitution frequency at the terminal position (row with pos == 1)."""
    with open(path) as fh:
        fh.readline()  # header: "pos\t5pC>T" (or similar)
        for line in fh:
            cols = line.split()
            if len(cols) >= 2:
                try:
                    return float(cols[1])
                except ValueError:
                    continue
    return 0.0


def authenticate(ct5: float, ga3: float, threshold: float) -> dict:
    present = ct5 >= threshold and ga3 >= threshold
    return {
        "ct_5prime_pos1": round(ct5, 4),
        "ga_3prime_pos1": round(ga3, 4),
        "damage_threshold": threshold,
        "damage_present": bool(present),
        "verdict": ("consistent with authentic ancient DNA"
                    if present else
                    "no clear terminal damage signal (modern contamination or low coverage?)"),
    }


def run(ct5_path: str, ga3_path: str, sample: str, threshold: float, out_json: str) -> dict:
    result = {"sample": sample,
              **authenticate(first_position_freq(ct5_path), first_position_freq(ga3_path),
                             threshold)}
    with open(out_json, "w") as fh:
        json.dump(result, fh, indent=2)
    return result


if "snakemake" in globals():  # pragma: no cover - exercised inside the workflow
    run(
        ct5_path=snakemake.input.ct5,           # noqa: F821
        ga3_path=snakemake.input.ga3,           # noqa: F821
        sample=snakemake.params.sample,         # noqa: F821
        threshold=float(snakemake.params.threshold),  # noqa: F821
        out_json=snakemake.output.json,         # noqa: F821
    )
