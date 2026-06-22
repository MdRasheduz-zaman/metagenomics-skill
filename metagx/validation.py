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

import glob
import gzip
import os
import re
import shutil
import subprocess
from typing import Dict, Iterable, List, Optional, Set, Tuple, Union

from . import formats


def blast_db_present(out_prefix: str) -> bool:
    """True if a BLAST nucleotide DB already exists at this prefix (idempotent build guard)."""
    return bool(glob.glob(out_prefix + "*.nin") or glob.glob(out_prefix + "*.nal"))


def normalize_seqid2taxid(map_path: str) -> Dict[str, str]:
    """kraken2 ``seqid2taxid.map`` -> ``{bare_accession: taxid}``.

    kraken2's map keys are the full internal seqid (``ACC|kraken:taxid|N``), but
    ``makeblastdb -taxid_map`` matches on the *bare* accession that ``-parse_seqids`` extracts —
    so the kraken2 map fails verbatim ("No sequences matched any of the taxids") and must be
    reduced to ``ACC -> taxid`` first. This is the bridge that makes the BLAST DB carry the
    EXACT taxids kraken2 assigned (synthetic or real), so validation is taxid-vs-taxid.
    """
    out: Dict[str, str] = {}
    with open(map_path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            acc = parts[0].split("|")[0].strip()
            taxid = parts[1].strip()
            if acc and taxid:
                out[acc] = taxid
    return out


def parse_names_dmp(path: str) -> Dict[str, str]:
    """taxid -> scientific name from an NCBI/kraken2 ``names.dmp`` (the in-sync name source —
    no NCBI taxdb needed, because both kraken2 and the taxid-tagged BLAST hits resolve here)."""
    # names.dmp columns: tax_id | name_txt | unique_name | name_class |
    out: Dict[str, str] = {}
    with open(path) as fh:
        for line in fh:
            f = [x.strip() for x in line.split("|")]
            if len(f) >= 4 and f[3] == "scientific name":
                out[f[0]] = f[1]
    return out


def kraken2_db_sources(db_dir: str) -> Dict[str, object]:
    """Resolve the inputs for an *in-sync* BLAST DB from a kraken2 DB directory.

    Returns the genome FASTA(s), the ``seqid2taxid.map``, and ``names.dmp`` when present.
    **Raises ValueError when the DB has no source genomes on disk** — a prebuilt/fetched index
    or a ``kraken2-build --clean``'d DB keeps only the opaque ``*.k2d`` hash, so there is
    nothing for blastn to align against; the caller must then supply the original genome FASTAs.
    """
    fastas: List[str] = []
    cl = os.path.join(db_dir, "custom_library.fasta")
    if os.path.isfile(cl):
        fastas = [cl]
    else:
        for pat in ("library/**/*.fna", "library/**/*.fa", "library/**/*.fasta",
                    "library/**/*.fna.gz", "library/**/*.fasta.gz"):
            fastas += glob.glob(os.path.join(db_dir, pat), recursive=True)
        fastas = sorted(set(fastas))
    if not fastas:
        raise ValueError(
            f"kraken2 DB '{db_dir}' has no source genomes on disk (a prebuilt/fetched index, or "
            f"a `kraken2-build --clean`'d build — only the *.k2d hash remains). blastn needs the "
            f"actual sequences: set validate.build_from to the genome FASTA(s) you built it from.")
    smap = os.path.join(db_dir, "seqid2taxid.map")
    names = os.path.join(db_dir, "taxonomy", "names.dmp")
    return {"fastas": fastas,
            "seqid2taxid": smap if os.path.isfile(smap) else None,
            "names_dmp": names if os.path.isfile(names) else None}


def build_blast_db(source: Union[str, List[str]], out_prefix: str, run: bool = True,
                   force: bool = False, taxid_map: Optional[str] = None) -> Dict[str, object]:
    """Build a BLAST+ nucleotide DB from a FASTA (or a folder of FASTAs) via makeblastdb.

    THIS is how the validation reference is kept *in scope* with the classifier: build it from
    the **same genomes** that built the kraken2/Bracken DB. Validating against a broader DB
    (e.g. full nt) measures a different benchmark — a read can match an organism the classifier
    never had a chance to call. ``-parse_seqids`` keeps subject titles so the organism is
    recoverable from ``stitle``; pass ``taxid_map`` (a bare-accession->taxid TSV, e.g. from
    ``normalize_seqid2taxid``) to tag subjects with kraken2's EXACT taxids so the comparison
    can be taxid-vs-taxid (no NCBI taxdb needed). Idempotent: skips if present.
    """
    result: Dict[str, object] = {"source": source, "db": os.path.abspath(out_prefix)}
    if not force and blast_db_present(out_prefix):
        result.update(ran=False, ok=True, skipped="already present")
        return result
    # a folder OR an explicit list of FASTAs -> concatenate into one input makeblastdb can read
    fasta = source
    files: List[str] = []
    if isinstance(source, (list, tuple)):
        files = [str(f) for f in source]
    elif os.path.isdir(source):
        files = [fp for fp in sorted(glob.glob(os.path.join(source, "*")))
                 if formats.read_format(fp) == "fasta"]
    if files:
        tmp_cat = out_prefix + ".sources.fasta"
        os.makedirs(os.path.dirname(out_prefix) or ".", exist_ok=True)
        with open(tmp_cat, "w") as out:
            for fp in files:
                opener = gzip.open if formats.is_gzipped(fp) else open
                with opener(fp, "rt") as fh:
                    shutil.copyfileobj(fh, out)
        fasta = tmp_cat
    cmd = ["makeblastdb", "-in", fasta, "-dbtype", "nucl", "-out", out_prefix, "-parse_seqids"]
    if taxid_map:
        cmd += ["-taxid_map", taxid_map]
    result["command"] = " ".join(cmd)
    if not run or not shutil.which("makeblastdb"):
        result["ran"] = False
        if not shutil.which("makeblastdb"):
            result["note"] = "makeblastdb not on PATH (install BLAST+, or use --use-conda)"
        return result
    os.makedirs(os.path.dirname(out_prefix) or ".", exist_ok=True)
    proc = subprocess.run(cmd, capture_output=True, text=True)
    result.update(ran=True, ok=(proc.returncode == 0 and blast_db_present(out_prefix)),
                  returncode=proc.returncode, taxid_mapped=bool(taxid_map),
                  tail=((proc.stdout or "") + (proc.stderr or ""))[-1000:])
    return result

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


def _blast_name(hit: Dict, name_resolver=None) -> str:
    """Best organism name for a BLAST hit, preferring the IN-SYNC taxonomy.

    When the BLAST DB was built with kraken2's taxid map, the hit carries kraken2's taxid in
    ``staxids`` and ``name_resolver`` (kraken2's names.dmp) turns it into the same name kraken2
    uses — an apples-to-apples comparison with no NCBI taxdb. Otherwise fall back to the BLAST
    scientific name, then the subject title (organism text in the FASTA header).
    """
    if name_resolver and hit.get("staxids"):
        tid = str(hit["staxids"]).split(";")[0].strip()
        nm = name_resolver(tid)
        if nm:
            return nm
    bname = (hit.get("sscinames") or "").strip()
    if not bname or bname.upper() == "N/A":
        bname = (hit.get("stitle") or "").strip()
    return bname


def assess(query_taxon: Dict[str, str], best_hits: Dict[str, Dict],
           level: str = "genus", name_resolver=None) -> Dict:
    """Per-query agreement between the classifier name and the best BLAST hit.

    Args:
      query_taxon  : {query_id: classifier_taxon_name} for the sequences we BLASTed.
      best_hits    : {query_id: blast6_row}.
      name_resolver: optional taxid->name (kraken2's names.dmp) — when the BLAST DB carries
                     kraken2's taxids, the comparison uses the SAME taxonomy (in-sync), not
                     BLAST's free-text title.
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
        bname = _blast_name(hit, name_resolver=name_resolver)
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
