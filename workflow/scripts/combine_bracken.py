"""Combine every per-(sample, sweep value) Bracken table into one long-format TSV.

Output columns: sample, label, name, taxonomy_id, taxonomy_lvl, new_est_reads,
fraction_total_reads. Filenames are {sample}.{label}.{level}.bracken.
"""

import os

import pandas as pd


def parse_name(path):
    base = os.path.basename(path)
    stem = base[: -len(".bracken")] if base.endswith(".bracken") else base
    sample, rest = stem.split(".", 1)
    label = rest.rsplit(".", 1)[0]  # drop trailing .<level>
    return sample, label


def main(tables, out_tsv):
    frames = []
    for path in tables:
        sample, label = parse_name(path)
        df = pd.read_csv(path, sep="\t")
        df.insert(0, "label", label)
        df.insert(0, "sample", sample)
        frames.append(df)
    combined = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    keep = [
        "sample", "label", "name", "taxonomy_id", "taxonomy_lvl",
        "new_est_reads", "fraction_total_reads",
    ]
    cols = [c for c in keep if c in combined.columns]
    combined = combined[cols] if cols else combined
    combined.to_csv(out_tsv, sep="\t", index=False)


if __name__ == "__main__":
    main(list(snakemake.input.tables), snakemake.output.tsv)  # noqa: F821
