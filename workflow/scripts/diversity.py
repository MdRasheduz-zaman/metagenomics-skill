"""Cross-sample diversity + ordination from the combined Bracken table.

Builds a sample x taxon matrix at one confidence label, then writes relative-abundance and
CLR matrices, an alpha-diversity table, a Bray-Curtis distance matrix, PCoA coordinates, a
top-taxa composition barplot, and (>=2 samples) a PCoA scatter.
"""

import csv
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from metagx import diversity as dv  # noqa: E402


def _write_tsv(path, header, rows):
    with open(path, "w") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(header)
        w.writerows(rows)


def main(bracken_combined, label, outdir, out_alpha, out_beta, out_pcoa, out_json,
         out_barplot, out_pcoa_png, top_n=15):
    os.makedirs(outdir, exist_ok=True)
    with open(bracken_combined) as fh:
        rows = [r for r in csv.DictReader(fh, delimiter="\t")
                if not label or r.get("label") in (label, None, "")]
    if not rows:  # fall back to all rows if the label filter emptied it
        with open(bracken_combined) as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))

    samples, taxa, mat = dv.build_matrix(rows)
    rel = dv.relative_abundance(mat)
    clr = dv.clr(mat)

    # relative-abundance + CLR matrices (taxa x samples)
    _write_tsv(os.path.join(outdir, "abundance_matrix.tsv"),
               ["taxon"] + samples,
               [[taxa[t]] + [round(rel[s, t], 6) for s in range(len(samples))]
                for t in range(len(taxa))])
    _write_tsv(os.path.join(outdir, "clr_matrix.tsv"),
               ["taxon"] + samples,
               [[taxa[t]] + [round(clr[s, t], 4) for s in range(len(samples))]
                for t in range(len(taxa))])

    alpha = dv.alpha_table(samples, mat)
    _write_tsv(out_alpha, ["sample", "richness", "shannon", "simpson", "pielou_evenness"],
               [[a["sample"], a["richness"], a["shannon"], a["simpson"], a["pielou_evenness"]]
                for a in alpha])

    bc = dv.braycurtis(mat)
    _write_tsv(out_beta, ["sample"] + samples,
               [[samples[i]] + [round(bc[i, j], 4) for j in range(len(samples))]
                for i in range(len(samples))])

    coords, explained = dv.pcoa(bc, n_axes=2)
    _write_tsv(out_pcoa, ["sample", "PCo1", "PCo2"],
               [[samples[i], round(coords[i, 0], 4), round(coords[i, 1], 4)]
                for i in range(len(samples))])

    # --- composition barplot (top taxa by mean relative abundance) ---
    if len(taxa):
        order = np.argsort(rel.mean(axis=0))[::-1][:top_n]
        plt.figure(figsize=(max(6, 1.2 * len(samples)), 6))
        bottom = np.zeros(len(samples))
        for t in order:
            plt.bar(samples, rel[:, t], bottom=bottom, label=taxa[t][:40])
            bottom += rel[:, t]
        plt.ylabel("relative abundance")
        plt.title("Community composition (top taxa)")
        plt.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=7)
        plt.xticks(rotation=45, ha="right")
        plt.tight_layout()
        plt.savefig(out_barplot, dpi=120)
        plt.close()

    # --- PCoA scatter (>=2 samples) ---
    if len(samples) >= 2:
        plt.figure(figsize=(6, 5))
        plt.scatter(coords[:, 0], coords[:, 1])
        for i, s in enumerate(samples):
            plt.annotate(s, (coords[i, 0], coords[i, 1]), fontsize=8)
        plt.xlabel(f"PCo1 ({explained[0]*100:.1f}%)")
        plt.ylabel(f"PCo2 ({explained[1]*100:.1f}%)")
        plt.title("PCoA (Bray-Curtis)")
        plt.tight_layout()
        plt.savefig(out_pcoa_png, dpi=120)
        plt.close()
    else:
        open(out_pcoa_png, "a").close()  # placeholder so the target exists

    with open(out_json, "w") as fh:
        json.dump({"n_samples": len(samples), "n_taxa": len(taxa), "label": label,
                   "alpha": alpha, "pcoa_explained": [round(float(x), 4) for x in explained]},
                  fh, indent=2)


if __name__ == "__main__":
    sm = snakemake  # noqa: F821
    main(sm.input.bracken, sm.params.label, sm.params.outdir,
         sm.output.alpha, sm.output.beta, sm.output.pcoa, sm.output.json,
         sm.output.barplot, sm.output.pcoa_png)
