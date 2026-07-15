"""Tag reference-FASTA headers with NCBI taxids for a ``taxonomy: real`` kraken2 build.

kraken2's real-taxonomy path reads each sequence's taxid from a ``>seqid|kraken:taxid|<taxid>``
header token; a custom reference usually lacks it, and ``metagx doctor`` only *warns*. This turns
that warning into a fix: rewrite the headers, resolving each sequence's accession -> NCBI taxid from
either a user-supplied TSV map (offline) or NCBI E-utilities (online), so the custom DB gets a real
lineage (species -> genus -> family -> ...), consistent with a standard kraken2 build.

Pure-Python, stdlib only (``urllib`` for the online path) — no new dependency, imports cleanly in
the CLI. The synthetic-taxonomy path (flat, offline, detection-only) stays as the escape hatch for
novel/air-gapped cases; this module serves the recommended real-taxonomy default.
"""
from __future__ import annotations

import json
import time
import urllib.parse
import urllib.request
from typing import Dict, List, Optional, Tuple

_ESUMMARY = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
_TAG = "kraken:taxid|"


def parse_accession(header: str) -> str:
    """Accession from a FASTA header: first whitespace token after ``>``, minus any pipe-delimited
    extras (including an existing ``|kraken:taxid|`` tag). ``>NC_001477.1 Dengue...`` ->
    ``NC_001477.1``; ``>NC_001477.1|kraken:taxid|11053 ...`` -> ``NC_001477.1``."""
    body = header.lstrip(">").strip()
    if not body:
        return ""
    return body.split()[0].split("|")[0]


def already_tagged(header: str) -> bool:
    return _TAG in header


def load_map(tsv_path: str) -> Dict[str, str]:
    """Read a 2-column ``accession<TAB>taxid`` map. Tolerant of a header row, ``#`` comments,
    blank lines, and extra columns. Keys are the accession *as written*; callers also try the
    version-stripped form."""
    out: Dict[str, str] = {}
    with open(tsv_path, encoding="utf-8-sig") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.replace("\t", " ").split()
            if len(parts) < 2:
                continue
            acc, taxid = parts[0], parts[1]
            if taxid.lower() in {"taxid", "tax_id"}:      # header row
                continue
            if taxid.isdigit():
                out[acc] = taxid
    return out


def _http_get(url: str, timeout: int = 30) -> str:
    req = urllib.request.Request(url, headers={"User-Agent": "metagx-tag-taxids"})
    with urllib.request.urlopen(req, timeout=timeout) as resp:      # noqa: S310 (fixed NCBI host)
        return resp.read().decode("utf-8", "replace")


def resolve_online(accessions: List[str], email: Optional[str] = None,
                   api_key: Optional[str] = None, batch: int = 180,
                   timeout: int = 30) -> Tuple[Dict[str, str], List[str]]:
    """Resolve accession -> taxid via NCBI E-utilities ``esummary`` (db=nuccore).

    Returns ``(mapping, failures)``. ``mapping`` is keyed by both the returned
    ``accessionversion`` and its version-stripped ``caption`` so either header form matches.
    Batched + rate-limited (NCBI: 3 req/s anonymous, 10 with an api_key). Network failures on a
    batch are collected in ``failures`` rather than raised, so a partial result is still usable."""
    mapping: Dict[str, str] = {}
    failures: List[str] = []
    uniq = list(dict.fromkeys(a for a in accessions if a))
    delay = 0.11 if api_key else 0.34
    for i in range(0, len(uniq), batch):
        chunk = uniq[i:i + batch]
        params = {"db": "nuccore", "id": ",".join(chunk), "retmode": "json"}
        if email:
            params["email"] = email
        if api_key:
            params["api_key"] = api_key
        url = f"{_ESUMMARY}?{urllib.parse.urlencode(params)}"
        try:
            data = json.loads(_http_get(url, timeout=timeout)).get("result", {})
        except Exception:                                  # noqa: BLE001 (network/JSON — keep partial)
            failures.extend(chunk)
            time.sleep(delay)
            continue
        for uid in data.get("uids", []):
            rec = data.get(uid, {})
            taxid = rec.get("taxid")
            if not taxid:
                continue
            accv = rec.get("accessionversion") or ""
            cap = rec.get("caption") or ""
            if accv:
                mapping[accv] = str(taxid)
            if cap:
                mapping.setdefault(cap, str(taxid))
        time.sleep(delay)
    return mapping, failures


def _lookup(acc: str, mapping: Dict[str, str]) -> Optional[str]:
    """Match an accession against the map, trying the exact then version-stripped form."""
    if acc in mapping:
        return mapping[acc]
    base = acc.split(".")[0]
    return mapping.get(base)


def tag_fasta(in_path: str, out_path: str, mapping: Dict[str, str],
              allow_missing: bool = False) -> Dict[str, object]:
    """Rewrite ``in_path`` headers to ``>acc|kraken:taxid|<taxid> <rest>`` using ``mapping``.

    Already-tagged headers are passed through untouched (idempotent). Returns a summary
    ``{n_seqs, n_tagged, n_already, n_missing, missing:[acc...]}``. Raises if any accession is
    unresolved and ``allow_missing`` is False — a silently-dropped sequence maps nowhere and
    quietly shrinks the DB, so this fails loud by default."""
    n_seqs = n_tagged = n_already = 0
    missing: List[str] = []
    with open(in_path, encoding="utf-8-sig") as fin, open(out_path, "w") as fout:
        for line in fin:
            if not line.startswith(">"):
                fout.write(line)
                continue
            n_seqs += 1
            if already_tagged(line):
                n_already += 1
                fout.write(line)
                continue
            acc = parse_accession(line)
            taxid = _lookup(acc, mapping)
            if taxid is None:
                missing.append(acc)
                fout.write(line)                            # leave untagged (reported below)
                continue
            rest = line[1:].rstrip("\n")
            first, _, tail = rest.partition(" ")
            fout.write(f">{first}|{_TAG}{taxid}" + (f" {tail}" if tail else "") + "\n")
            n_tagged += 1
    summary = {"n_seqs": n_seqs, "n_tagged": n_tagged, "n_already": n_already,
               "n_missing": len(missing), "missing": missing[:50]}
    if missing and not allow_missing:
        raise ValueError(
            f"{len(missing)} sequence(s) had no taxid (e.g. {', '.join(missing[:5])}). Provide them "
            f"in --map, use --online to resolve, or --allow-missing to skip them (they will NOT be "
            f"classifiable).")
    return summary
