"""Reconcile contig-level taxonomy with read-level classification.

Joins three sources for one sample:
  * kraken2 per-contig calls   (taxonomy of each assembled contig, --use-names)
  * per-contig depth           (coverage, from jgi_summarize_bam_contig_depths)
  * read-level kraken2 report  (what the raw reads said, at one confidence)

Produces:
  * <sample>.contig_taxonomy.tsv  — per contig: length, coverage, taxid, taxon, C/U
  * <sample>.reconciliation.tsv   — per taxon: read abundance vs contig evidence + concordance
  * <sample>.flags.tsv            — discordances worth a look (contig-only taxa, high-cov
                                    unclassified contigs = candidate novel/divergent/chimeric)
  * <sample>.reconcile.json       — compact summary for the report / LLM

Interpretation:
  read abundance answers "how much"; contig evidence answers "what is confidently present
  and reconstructable". coverage-weighted contig abundance (length x depth) is the bridge:
  taxonomy from the long consensus, quantity from the reads. Taxa are joined by taxid;
  note contig LCA calls can sit at a higher rank than the read/species level.
"""

import csv
import json
import re

_TAXID_RE = re.compile(r"taxid\s+(\d+)")


def parse_contig_calls(path):
    """kraken2 --use-names per-sequence output: C/U, id, 'name (taxid N)', length, lca."""
    rows = []
    with open(path) as fh:
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) < 4:
                continue
            status, contig, name_field, length = c[0], c[1], c[2], c[3]
            m = _TAXID_RE.search(name_field)
            taxid = int(m.group(1)) if m else 0
            name = name_field.split(" (taxid")[0].strip() if "(taxid" in name_field else name_field
            try:
                length = int(str(length).split("|")[0])  # paired lengths "a|b"; contigs are single
            except ValueError:
                length = 0
            rows.append({"contig": contig, "status": status, "taxid": taxid,
                         "taxon": name, "length": length})
    return rows


def parse_depth(path):
    """jgi depth file: contigName, contigLen, totalAvgDepth, <per-bam>..."""
    depth = {}
    try:
        with open(path) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                try:
                    depth[r["contigName"]] = float(r.get("totalAvgDepth", 0) or 0)
                except (ValueError, KeyError):
                    depth[r["contigName"]] = 0.0
    except OSError:
        pass
    return depth


def parse_cat(path):
    """CAT add_names output -> {contig: {taxid, taxon}} (per-ORF voting classification)."""
    out = {}
    if not path or not __import__("os").path.exists(path):
        return out
    with open(path) as fh:
        fh.readline()  # header
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) < 4:
                continue
            taxid = 0
            for tok in reversed(c[3].replace("*", "").split(";")):
                tok = tok.strip()
                if tok.isdigit():
                    taxid = int(tok)
                    break
            name = ""
            for col in reversed(c[5:]):
                col = col.strip()
                if col and col != "NA":
                    name = col.split(":")[0].strip()
                    break
            out[c[0]] = {"taxid": taxid, "taxon": name}
    return out


def parse_breadth(path):
    """samtools coverage output -> {contig: percent_bases_covered}."""
    out = {}
    try:
        with open(path) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                key = r.get("#rname") or r.get("rname")
                try:
                    out[key] = float(r.get("coverage", 0) or 0)  # 'coverage' = % bases covered
                except ValueError:
                    out[key] = 0.0
    except OSError:
        pass
    return out


def parse_read_report(path):
    """kraken2 read report -> {taxid: {name, pct, reads}} for species (rank S) rows."""
    reads = {}
    with open(path) as fh:
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) < 6 or c[3].strip() != "S":
                continue
            try:
                taxid = int(c[4])
            except ValueError:
                continue
            reads[taxid] = {"name": c[5].strip(),
                            "pct": float(c[0]), "reads": int(c[1])}
    return reads


def main(contig_calls, depth_path, read_report, sample, read_label,
         out_contig, out_taxon, out_flags, out_json, cat_named="", breadth_path=""):
    contigs = parse_contig_calls(contig_calls)
    depth = parse_depth(depth_path)
    breadth = parse_breadth(breadth_path)
    reads = parse_read_report(read_report)
    cat = parse_cat(cat_named)  # {} if CAT not run

    # ---- per-contig table (kraken2 vs CAT methods) ----
    n_cat_classified = n_agree = n_conflict = 0
    for row in contigs:
        row["coverage"] = round(depth.get(row["contig"], 0.0), 3)
        row["breadth_pct"] = round(breadth.get(row["contig"], 0.0), 2)
        ct = cat.get(row["contig"], {})
        row["cat_taxid"] = ct.get("taxid", 0)
        row["cat_taxon"] = ct.get("taxon", "")
        k_ok = row["status"] == "C" and row["taxid"]
        c_ok = bool(row["cat_taxid"])
        if c_ok:
            n_cat_classified += 1
        if k_ok and c_ok:
            if row["taxid"] == row["cat_taxid"]:
                row["methods_agree"] = "yes"; n_agree += 1
            else:
                row["methods_agree"] = "no"; n_conflict += 1
        else:
            row["methods_agree"] = "partial" if (k_ok or c_ok) else ""
    with open(out_contig, "w") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["contig", "length", "coverage", "breadth_pct", "taxid", "taxon",
                    "classified", "cat_taxid", "cat_taxon", "methods_agree"])
        for r in sorted(contigs, key=lambda x: x["coverage"], reverse=True):
            w.writerow([r["contig"], r["length"], r["coverage"], r["breadth_pct"],
                        r["taxid"], r["taxon"], r["status"], r["cat_taxid"],
                        r["cat_taxon"], r["methods_agree"]])

    # ---- aggregate contig side by taxid (classified only) ----
    contig_tax = {}
    total_covwt = 0.0
    for r in contigs:
        if r["status"] != "C" or r["taxid"] == 0:
            continue
        covwt = r["length"] * max(r["coverage"], 0.0)
        total_covwt += covwt
        t = contig_tax.setdefault(r["taxid"], {"name": r["taxon"], "n_contigs": 0,
                                               "assembled_bp": 0, "covwt": 0.0})
        t["n_contigs"] += 1
        t["assembled_bp"] += r["length"]
        t["covwt"] += covwt

    # ---- per-taxon reconciliation (union of read + contig taxids) ----
    all_taxids = set(reads) | set(contig_tax)
    taxon_rows = []
    for taxid in all_taxids:
        rd = reads.get(taxid)
        ct = contig_tax.get(taxid)
        in_reads, in_contigs = rd is not None, ct is not None
        concordance = "both" if (in_reads and in_contigs) else ("reads_only" if in_reads else "contigs_only")
        taxon_rows.append({
            "taxid": taxid,
            "taxon": (rd or ct)["name"],
            "read_pct": round(rd["pct"], 4) if rd else 0.0,
            "read_reads": rd["reads"] if rd else 0,
            "n_contigs": ct["n_contigs"] if ct else 0,
            "assembled_bp": ct["assembled_bp"] if ct else 0,
            "cov_weighted_pct": round(100.0 * ct["covwt"] / total_covwt, 4) if (ct and total_covwt) else 0.0,
            "concordance": concordance,
        })
    taxon_rows.sort(key=lambda x: (x["read_pct"], x["cov_weighted_pct"]), reverse=True)
    with open(out_taxon, "w") as fh:
        cols = ["taxid", "taxon", "read_pct", "read_reads", "n_contigs",
                "assembled_bp", "cov_weighted_pct", "concordance"]
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(taxon_rows)

    # ---- discordance / chimera flags ----
    flags = []
    for r in taxon_rows:
        if r["concordance"] == "contigs_only":
            flags.append({"type": "contig_only_taxon", "taxid": r["taxid"], "taxon": r["taxon"],
                          "note": "assembled & classified but not seen in reads at species level "
                                  "(novel-to-reads, higher-rank read call, or possible misassembly)"})
        elif r["concordance"] == "reads_only" and r["read_pct"] >= 1.0:
            flags.append({"type": "reads_only_taxon", "taxid": r["taxid"], "taxon": r["taxon"],
                          "note": f"{r['read_pct']}% of reads but did not assemble "
                                  "(low coverage / fragmented / strain mixture)"})
    # high-coverage unclassified contigs = candidate divergent/novel or chimeric
    unc = sorted((r for r in contigs if r["status"] != "C" or r["taxid"] == 0),
                 key=lambda x: x["coverage"], reverse=True)[:10]
    for r in unc:
        if r["coverage"] > 0:
            flags.append({"type": "unclassified_contig", "contig": r["contig"],
                          "length": r["length"], "coverage": r["coverage"],
                          "note": "unclassified contig with coverage — candidate novel/divergent or chimeric"})
    with open(out_flags, "w") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["type", "taxid_or_contig", "taxon_or_length", "coverage_or_pct", "note"])
        for f in flags:
            key = f.get("taxid", f.get("contig", ""))
            val = f.get("taxon", f.get("length", ""))
            extra = f.get("coverage", f.get("read_pct", ""))
            w.writerow([f["type"], key, val, extra, f["note"]])

    # ---- compact summary ----
    counts = {"both": 0, "reads_only": 0, "contigs_only": 0}
    for r in taxon_rows:
        counts[r["concordance"]] += 1
    summary = {
        "sample": sample,
        "read_label": read_label,
        "n_contigs": len(contigs),
        "n_contigs_classified": sum(1 for r in contigs if r["status"] == "C" and r["taxid"]),
        "assembled_bp": sum(r["length"] for r in contigs),
        "taxa_concordance": counts,
        "n_flags": len(flags),
        "top_taxa": taxon_rows[:10],
    }
    if cat:
        summary["cat_cross_check"] = {
            "n_cat_classified": n_cat_classified,
            "kraken2_cat_agree": n_agree,
            "kraken2_cat_conflict": n_conflict,
        }
    with open(out_json, "w") as fh:
        json.dump(summary, fh, indent=2)


if __name__ == "__main__":
    sm = snakemake  # noqa: F821
    main(sm.input.contig_calls, sm.input.depth, sm.input.read_report,
         sm.params.sample, sm.params.read_label,
         sm.output.per_contig, sm.output.per_taxon, sm.output.flags, sm.output.summary,
         cat_named=sm.params.cat_named, breadth_path=sm.input.breadth)
