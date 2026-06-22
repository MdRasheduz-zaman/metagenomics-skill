"""BLAST-based validation of taxonomic classifications (the `validate` module's core logic).

A kraken2/Bracken assignment is a k-mer call, not proof. This module takes the classifier's
top taxa, pulls a representative set of sequences for each, BLASTs them against a reference,
and asks the only question that matters for validation: *does the best alignment's organism
agree with what the classifier said?* It reports a per-taxon agreement rate and an overall
concordance, so a user can see which calls are corroborated by alignment and which are likely
DB-completeness artifacts.

Pure-Python and dependency-light (stdlib only): parsing + the agreement check live here and
are unit-tested without BLAST installed; the workflow rule does the IO (read extraction) and
shells out to `blastn`. Matches the project's no-scipy/no-R convention.
"""
from __future__ import annotations

import gzip
import re
from typing import Dict, Iterable, List, Optional, Set, Tuple

from . import formats

# The outfmt-6 columns the validate rule requests. Order matters — it's how we parse.
BLAST6_FIELDS = ["qseqid", "sseqid", "pident", "length", "evalue", "bitscore",
                 "staxids", "sscinames", "stitle"]

_TAXID_IN_NAME = re.compile(r"\(taxid\s+(\d+)\)")
_NONWORD = re.compile(r"[^A-Za-z0-9]+")


def parse_kraken_assignments(text: str) -> Dict[str, str]:
    """Map read id -> assigned taxid from a kraken2 per-read `.kraken` output.

    Each line is: ``C|U <readid> <taxid-or-"name (taxid N)"> <len> <lca>``. Handles both the
    bare-taxid form (kraken2 default) and the --use-names form. Unclassified (U / taxid 0) reads
    are skipped — there's nothing to validate.
    """
    out: Dict[str, str] = {}
    for line in text.splitlines():
        parts = line.split("\t")
        if len(parts) < 3 or parts[0] != "C":
            continue
        readid, taxfield = parts[1], parts[2]
        m = _TAXID_IN_NAME.search(taxfield)
        taxid = m.group(1) if m else taxfield.strip()
        if taxid and taxid != "0":
            out[readid] = taxid
    return out


def top_taxa(kreport_path: str, level: str = "S", top_n: int = 10) -> List[Dict]:
    """Top taxa at a rank from a kraken2 report, by reads-rooted-here (kreport column 2).

    Returns [{taxid, name, rank, reads}] sorted descending. Uses ``formats.kreport_row`` so it
    survives the --report-minimizer-data column shift (never fixed indices)."""
    want = level.upper()
    rows: List[Dict] = []
    with open(kreport_path) as fh:
        for line in fh:
            row = formats.kreport_row(line)
            if not row or not row.get("taxid"):
                continue
            rank = str(row.get("rank", "")).upper()
            # accept exact rank or a sub-rank starting with it (S, S1, S2 ...) when level is S/G
            if rank != want and not (len(want) == 1 and rank.startswith(want)):
                continue
            try:
                reads = int(row["clade_reads"])  # reads rooted at this clade
            except (ValueError, KeyError, TypeError):
                reads = 0
            rows.append({"taxid": str(row["taxid"]), "name": str(row.get("name", "")).strip(),
                         "rank": rank, "reads": reads})
    rows.sort(key=lambda r: r["reads"], reverse=True)
    return rows[:top_n]


def _norm_id(raw: str) -> str:
    """A read id as kraken records it: drop a leading @/>, trim at whitespace, drop /1 /2."""
    rid = raw.lstrip("@>").split()[0] if raw.strip() else ""
    if rid.endswith(("/1", "/2")):
        rid = rid[:-2]
    return rid


def extract_sequences(paths: Iterable[str], wanted: Set[str]) -> Dict[str, str]:
    """Return {read_id: sequence} for ids in ``wanted`` from FASTA/FASTQ(.gz) files.

    Pure-Python, streams the files, stops once every wanted id is found. Read ids are
    normalized (`_norm_id`) so paired mate suffixes and description fields don't prevent a
    match against kraken's per-read ids.
    """
    found: Dict[str, str] = {}
    remaining = set(wanted)
    for path in paths:
        if not remaining:
            break
        opener = gzip.open if formats.is_gzipped(str(path)) else open
        is_fa = formats.read_format(str(path)) == "fasta"
        with opener(path, "rt") as fh:
            if is_fa:
                rid, seq = None, []
                for line in fh:
                    if line.startswith(">"):
                        if rid is not None and rid in remaining:
                            found[rid] = "".join(seq); remaining.discard(rid)
                        rid, seq = _norm_id(line), []
                    else:
                        seq.append(line.strip())
                if rid is not None and rid in remaining:
                    found[rid] = "".join(seq); remaining.discard(rid)
            else:  # fastq: 4-line records
                while True:
                    header = fh.readline()
                    if not header:
                        break
                    seq = fh.readline().strip()
                    fh.readline(); fh.readline()  # '+' and qual
                    rid = _norm_id(header)
                    if rid in remaining:
                        found[rid] = seq; remaining.discard(rid)
                        if not remaining:
                            break
    return found


def write_fasta(seqs: Dict[str, str], path: str) -> int:
    """Write {id: seq} to a FASTA file; return the count written."""
    n = 0
    with open(path, "w") as fh:
        for rid, seq in seqs.items():
            if seq:
                fh.write(f">{rid}\n{seq}\n"); n += 1
    return n


def parse_blast6(text: str, fields: Optional[List[str]] = None) -> List[Dict]:
    """Parse BLAST outfmt-6 (tab-separated) into dicts keyed by ``fields``."""
    fields = fields or BLAST6_FIELDS
    rows: List[Dict] = []
    for line in text.splitlines():
        if not line.strip():
            continue
        vals = line.split("\t")
        row = dict(zip(fields, vals))
        for k in ("pident", "evalue", "bitscore"):
            if k in row:
                try:
                    row[k] = float(row[k])
                except ValueError:
                    row[k] = 0.0
        rows.append(row)
    return rows


def best_hit_per_query(rows: List[Dict]) -> Dict[str, Dict]:
    """Highest-bitscore hit per query id."""
    best: Dict[str, Dict] = {}
    for r in rows:
        q = r.get("qseqid")
        if q is None:
            continue
        if q not in best or r.get("bitscore", 0.0) > best[q].get("bitscore", 0.0):
            best[q] = r
    return best


def _tokens(name: str) -> List[str]:
    return [t for t in _NONWORD.split(name.lower()) if t]


def names_agree(classifier_name: str, blast_name: str, level: str = "genus") -> bool:
    """Does a BLAST subject organism corroborate the classifier's taxon name?

    genus level => share the first (genus) token; species level => share the first two
    (genus + species epithet). Robust to the classifier's "Genus species" vs BLAST's free-text
    stitle/sscinames by token-overlap rather than exact string equality.
    """
    c, b = _tokens(classifier_name), _tokens(blast_name)
    if not c or not b:
        return False
    if level == "species":
        return len(c) >= 2 and c[0] == b[0] and c[1] in b
    return c[0] == b[0] or c[0] in b  # genus token present in the BLAST hit


def assess(query_taxon: Dict[str, str], best_hits: Dict[str, Dict],
           level: str = "genus") -> Dict:
    """Per-query agreement between the classifier name and the best BLAST hit.

    Args:
      query_taxon : {query_id: classifier_taxon_name} for the sequences we BLASTed.
      best_hits   : {query_id: blast6_row} (use sscinames, falling back to stitle).
    Returns aggregate counts + per-query verdicts.
    """
    n = len(query_taxon)
    with_hits = 0
    agree = 0
    per_query: List[Dict] = []
    for q, cname in query_taxon.items():
        hit = best_hits.get(q)
        if not hit:
            per_query.append({"query": q, "classifier": cname, "blast": None,
                              "has_hit": False, "agree": False})
            continue
        with_hits += 1
        # Prefer the scientific name; but without NCBI's taxdb installed blastn writes "N/A"
        # there, so fall back to the subject title (which carries the organism in its text).
        bname = (hit.get("sscinames") or "").strip()
        if not bname or bname.upper() == "N/A":
            bname = (hit.get("stitle") or "").strip()
        ok = names_agree(cname, bname, level=level)
        agree += int(ok)
        per_query.append({"query": q, "classifier": cname, "blast": bname,
                          "pident": hit.get("pident"), "has_hit": True, "agree": ok})
    return {
        "n_queries": n,
        "n_with_hits": with_hits,
        "n_agree": agree,
        "hit_rate": round(with_hits / n, 4) if n else 0.0,
        # agreement is over queries that produced a hit (no-hit ≠ disagreement)
        "agreement_rate": round(agree / with_hits, 4) if with_hits else 0.0,
        "per_query": per_query,
    }


def verdict(agreement_rate: float, hit_rate: float) -> str:
    """A one-word call the report/advisor can surface."""
    if hit_rate < 0.2:
        return "inconclusive"  # too few alignments to judge (DB too small / too divergent)
    if agreement_rate >= 0.8:
        return "corroborated"
    if agreement_rate >= 0.5:
        return "partial"
    return "discordant"


def summarize(sample: str, per_taxon: List[Dict], level: str, db: str) -> Dict:
    """Roll per-taxon assessments into the sample's validation report."""
    tot_q = sum(t["n_queries"] for t in per_taxon)
    tot_hit = sum(t["n_with_hits"] for t in per_taxon)
    tot_agree = sum(t["n_agree"] for t in per_taxon)
    overall_agree = round(tot_agree / tot_hit, 4) if tot_hit else 0.0
    overall_hit = round(tot_hit / tot_q, 4) if tot_q else 0.0
    return {
        "sample": sample,
        "level": level,
        "db": db,
        "n_taxa_validated": len(per_taxon),
        "n_queries": tot_q,
        "n_with_hits": tot_hit,
        "n_agree": tot_agree,
        "overall_hit_rate": overall_hit,
        "overall_agreement_rate": overall_agree,
        "verdict": verdict(overall_agree, overall_hit),
        "taxa": per_taxon,
    }
