"""Differential abundance between two sample groups, from the combined Bracken table.

Builds a sample x taxon matrix at one confidence label, maps samples to their `group` label
(sample-sheet column), and runs a CLR + permutation test (metagx.differential). Writes a
ranked results table, a JSON summary, and a volcano plot. Pure-Python (numpy/matplotlib).
"""

import csv
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from metagx import diversity as dv  # noqa: E402
from metagx import differential as da  # noqa: E402


def main(bracken_combined, label, groups, out_tsv, out_json, out_png,
         n_perm=999, fdr=0.05, seed=42, reference_group=None):
    os.makedirs(os.path.dirname(out_tsv), exist_ok=True)
    with open(bracken_combined) as fh:
        rows = [r for r in csv.DictReader(fh, delimiter="\t")
                if not label or r.get("label") in (label, None, "")]
    if not rows:
        with open(bracken_combined) as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))

    samples, taxa, mat = dv.build_matrix(rows)
    res_rows, summary = da.differential_abundance(
        samples, taxa, mat, groups, n_perm=n_perm, fdr=fdr, seed=seed,
        reference_group=reference_group)

    header = list(res_rows[0].keys()) if res_rows else [
        "taxon", "clr_diff", "effect_size", "p_value", "q_value", "significant"]
    with open(out_tsv, "w") as fh:
        w = csv.DictWriter(fh, fieldnames=header, delimiter="\t")
        w.writeheader()
        w.writerows(res_rows)

    with open(out_json, "w") as fh:
        json.dump({"summary": summary,
                   "significant_taxa": [r["taxon"] for r in res_rows if r["significant"]]},
                  fh, indent=2)

    # --- volcano: effect size (x) vs -log10(q) (y), significant points highlighted ---
    if res_rows:
        x = np.array([r["clr_diff"] for r in res_rows])
        q = np.array([max(r["q_value"], 1e-6) for r in res_rows])
        y = -np.log10(q)
        sig = np.array([r["significant"] for r in res_rows])
        plt.figure(figsize=(6, 5))
        plt.scatter(x[~sig], y[~sig], s=12, c="#bbbbbb", label="ns")
        if sig.any():
            plt.scatter(x[sig], y[sig], s=18, c="#cc3333", label=f"q<{summary['fdr']}")
        plt.axhline(-np.log10(summary["fdr"]), ls="--", c="#888888", lw=0.8)
        plt.xlabel(f"CLR difference ({summary['group_a']} − {summary['group_b']})")
        plt.ylabel("−log10(q)")
        plt.title("Differential abundance (volcano)")
        plt.legend(fontsize=8)
        plt.tight_layout()
        plt.savefig(out_png, dpi=120)
        plt.close()
    else:
        open(out_png, "a").close()


if __name__ == "__main__":
    sm = snakemake  # noqa: F821
    d = sm.params.differential
    main(sm.input.bracken, sm.params.label, sm.params.groups,
         sm.output.tsv, sm.output.json, sm.output.png,
         n_perm=d.get("n_permutations", 999), fdr=d.get("fdr", 0.05),
         seed=d.get("seed", 42), reference_group=d.get("reference_group"))
