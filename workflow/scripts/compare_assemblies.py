"""Compare unfiltered vs taxonomically-filtered assemblies for one sample.

Reports contig count, total bp, N50, and longest contig for each, plus deltas — so you can
SEE whether filtering improved recovery (e.g. host depletion → longer contigs) or just
discarded data (fewer/shorter contigs = filtering hurt the target).
"""

import json


def stats(fa):
    lens, cur = [], 0
    try:
        with open(fa) as fh:
            for line in fh:
                if line.startswith(">"):
                    if cur:
                        lens.append(cur)
                    cur = 0
                else:
                    cur += len(line.strip())
        if cur:
            lens.append(cur)
    except OSError:
        pass
    lens.sort(reverse=True)
    total = sum(lens)
    n50, c = 0, 0
    for L in lens:
        c += L
        if c >= total / 2:
            n50 = L
            break
    return {"contigs": len(lens), "total_bp": total,
            "longest": lens[0] if lens else 0, "n50": n50}


def main(unfiltered, filtered, sample, out_tsv, out_json):
    u, f = stats(unfiltered), stats(filtered)
    metrics = ["contigs", "total_bp", "longest", "n50"]
    with open(out_tsv, "w") as fh:
        fh.write("metric\tunfiltered\tfiltered\tdelta\n")
        for m in metrics:
            fh.write(f"{m}\t{u[m]}\t{f[m]}\t{f[m] - u[m]}\n")
    with open(out_json, "w") as fh:
        json.dump({"sample": sample, "unfiltered": u, "filtered": f,
                   "delta": {m: f[m] - u[m] for m in metrics}}, fh, indent=2)


if __name__ == "__main__":
    sm = snakemake  # noqa: F821
    main(sm.input.unfiltered, sm.input.filtered, sm.params.sample,
         sm.output.tsv, sm.output.json)
