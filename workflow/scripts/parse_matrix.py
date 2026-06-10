"""Build a per-sample comparison matrix across the sweep axis.

Reads the kraken2 reports for one sample at every sweep value, extracts species-level
read counts, and writes (1) a compact JSON matrix for the LLM and (2) a line plot
showing how each top organism's assigned reads change across the sweep.

Robust to the extra columns added by kraken2 --report-minimizer-data: rank code,
taxid, and name are always read from the last three columns.
"""

import json
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402


def species_counts(report_path):
    df = pd.read_csv(report_path, sep="\t", header=None, dtype=str)
    # reads-in-clade is column 1; rank/taxid/name are the last three columns.
    reads = pd.to_numeric(df[1], errors="coerce").fillna(0).astype(int)
    rank = df.iloc[:, -3].str.strip()
    name = df.iloc[:, -1].str.strip()
    out = pd.DataFrame({"name": name, "reads": reads})[rank == "S"]
    return out.set_index("name")["reads"]


def main(reports, labels, sweep_param, sample, out_json, out_png, top_n=15):
    matrix = pd.DataFrame()
    for path, label in zip(reports, labels):
        matrix[label] = species_counts(path)
    matrix = matrix.fillna(0).astype(int)

    matrix["__total__"] = matrix.sum(axis=1)
    top = matrix.sort_values("__total__", ascending=False).head(top_n).drop(columns="__total__")

    payload = {
        "sample": sample,
        "sweep_param": sweep_param,
        "sweep_labels": list(labels),
        "n_species_total": int((matrix.drop(columns="__total__") > 0).any(axis=1).sum()),
        "top_species": top.to_dict(orient="index"),
    }
    with open(out_json, "w") as fh:
        json.dump(payload, fh, indent=2)

    plt.figure(figsize=(11, 6))
    for organism in top.index:
        plt.plot(list(top.columns), top.loc[organism], marker="o", label=organism)
    plt.title(f"{sample}: classification sensitivity across {sweep_param}")
    plt.ylabel("Reads assigned (clade)")
    plt.xlabel(sweep_param)
    plt.xticks(rotation=45, ha="right")
    plt.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=8)
    plt.tight_layout()
    plt.savefig(out_png, dpi=120)
    plt.close()


if __name__ == "__main__":
    main(
        list(snakemake.input.reports),  # noqa: F821
        list(snakemake.params.labels),  # noqa: F821
        snakemake.params.sweep_param,  # noqa: F821
        snakemake.params.sample,  # noqa: F821
        snakemake.output.json,  # noqa: F821
        snakemake.output.png,  # noqa: F821
    )
