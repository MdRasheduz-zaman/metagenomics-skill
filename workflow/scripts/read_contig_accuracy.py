"""Score read-level kraken2 calls against the (more reliable) contig-level call.

For every read that maps to a contig, compare what kraken2 said about the *read* to what it
said about the *contig the read assembled into*. The contig call (longer, error-corrected
consensus) is used as a practical ground truth. Each aligned read is bucketed:

  concordant_exact    read taxid == contig taxid
  concordant_lineage  one is an ancestor of the other (agree at a coarser/finer rank)
  discordant          read classified to a different lineage (effectively misclassified)
  read_unclassified   contig classified, read was not
  contig_unclassified contig itself unclassified (can't judge the read) — excluded from rates

Inputs: the sorted BAM (reads→contigs), the read-level .kraken (col3 = taxid), the
contig-level .kraken (--use-names, col3 = "name (taxid N)"), and nodes.dmp for lineage.
No new tools — just samtools (already used for mapping) + stdlib.
"""

import csv
import json
import re
import subprocess
from collections import Counter, defaultdict

_TAXID_RE = re.compile(r"taxid\s+(\d+)")
_CATS = ["concordant_exact", "concordant_lineage", "discordant",
         "read_unclassified", "contig_unclassified"]


def load_nodes(path):
    parent = {}
    try:
        with open(path) as fh:
            for line in fh:
                p = [x.strip() for x in line.split("|")]
                if len(p) >= 2 and p[0].isdigit():
                    parent[int(p[0])] = int(p[1])
    except OSError:
        pass
    return parent


def _is_ancestor(anc, desc, parent):
    x, seen = desc, 0
    while x in parent and parent[x] != x and seen < 200:
        x = parent[x]
        if x == anc:
            return True
        seen += 1
    return False


def _related(a, b, parent):
    return _is_ancestor(a, b, parent) or _is_ancestor(b, a, parent)


def parse_contig_taxa(path):
    out = {}
    with open(path) as fh:
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) < 3:
                continue
            m = _TAXID_RE.search(c[2])
            out[c[1]] = int(m.group(1)) if m else 0
    return out


def parse_read_taxa(path):
    out = {}
    with open(path) as fh:
        for line in fh:
            c = line.split("\t")
            if len(c) < 3:
                continue
            try:
                out[c[1]] = int(c[2])
            except ValueError:
                out[c[1]] = 0
    return out


def main(bam, read_kraken, contig_kraken, nodes_path, sample, out_tsv, out_json):
    contig_tax = parse_contig_taxa(contig_kraken)
    read_tax = parse_read_taxa(read_kraken)
    parent = load_nodes(nodes_path)

    per_contig = defaultdict(Counter)
    overall = Counter()
    # primary, mapped alignments only (-F 0x904: unmapped, secondary, supplementary)
    proc = subprocess.Popen(["samtools", "view", "-F", "0x904", bam],
                            stdout=subprocess.PIPE, text=True)
    for line in proc.stdout:
        f = line.split("\t", 3)
        if len(f) < 3:
            continue
        qname, rname = f[0], f[2]
        ctax = contig_tax.get(rname, 0)
        rtax = read_tax.get(qname, 0)
        if ctax == 0:
            cat = "contig_unclassified"
        elif rtax == 0:
            cat = "read_unclassified"
        elif rtax == ctax:
            cat = "concordant_exact"
        elif _related(rtax, ctax, parent):
            cat = "concordant_lineage"
        else:
            cat = "discordant"
        per_contig[rname][cat] += 1
        overall[cat] += 1
    proc.wait()

    with open(out_tsv, "w") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["contig", "contig_taxid", "total_reads"] + _CATS)
        for contig in sorted(per_contig, key=lambda c: sum(per_contig[c].values()), reverse=True):
            cc = per_contig[contig]
            w.writerow([contig, contig_tax.get(contig, 0), sum(cc.values())]
                       + [cc.get(k, 0) for k in _CATS])

    judged = sum(overall[k] for k in _CATS if k != "contig_unclassified")
    concordant = overall["concordant_exact"] + overall["concordant_lineage"]
    summary = {
        "sample": sample,
        "reads_aligned": sum(overall.values()),
        "reads_on_classified_contigs": judged,
        "counts": {k: overall.get(k, 0) for k in _CATS},
        "rates_on_classified_contigs": {
            "concordant_pct": round(100.0 * concordant / judged, 2) if judged else 0.0,
            "discordant_pct": round(100.0 * overall["discordant"] / judged, 2) if judged else 0.0,
            "read_unclassified_pct": round(100.0 * overall["read_unclassified"] / judged, 2) if judged else 0.0,
        },
    }
    with open(out_json, "w") as fh:
        json.dump(summary, fh, indent=2)


if __name__ == "__main__":
    sm = snakemake  # noqa: F821
    main(sm.input.bam, sm.input.read_kraken, sm.input.contig_kraken, sm.params.nodes,
         sm.params.sample, sm.output.tsv, sm.output.json)
