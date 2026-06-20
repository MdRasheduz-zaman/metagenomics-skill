"""`metagx probe` — measured pre-flight context from the user's own reads.

Profiles a bounded head subsample of *every* sample in the sheet, locally, and emits only
non-reconstructive aggregate statistics (length / quality / GC / duplication / inferred
platform). The output `context` dict feeds ``registry.interview_spec(context=)`` so
goal/data-conditional promotion fires on measured facts, not asserted ones.

Constraints (see docs/spec-probe.md): consent-gated, local-only, never stores read sequences
or IDs, pure-Python (no scipy, no bio tools on PATH). The one reference-dependent metric in the
spec — host fraction — is out of scope for this MVP and reported as null.
"""

from __future__ import annotations

import csv
import gzip
import os
from typing import Any, Dict, List, Optional

from . import consent
from .formats import is_gzipped, read_format

# Inference thresholds (deliberately conservative; tune from evidence/*.yaml later).
LONG_MIN_LEN = 400      # median length above this => long-read
HIFI_MAX_ERR = 0.02     # long + estimated error below this => accurate long-read class
LOW_Q20_FRAC = 0.90     # below this fraction of Q>=20 bases => flag low quality
GZIP_RATIO = 4.0        # rough fastq compression ratio for size extrapolation
PREFIX_BP = 32          # duplication probe: hash of the first PREFIX_BP bases (never stored raw)

_SHORT = {"illumina", "mgi", "bgi"}


def _open(path: str):
    return gzip.open(path, "rt") if is_gzipped(path) else open(path, "rt")


def _median(xs: List[float]) -> float:
    if not xs:
        return 0
    s = sorted(xs)
    n = len(s)
    return s[n // 2] if n % 2 else (s[n // 2 - 1] + s[n // 2]) / 2


def _percentile(xs: List[int], p: float) -> int:
    if not xs:
        return 0
    s = sorted(xs)
    return s[min(len(s) - 1, int(p * (len(s) - 1)))]


# --------------------------------------------------------------------------- #
# Sample sheet                                                                 #
# --------------------------------------------------------------------------- #
def load_sheet(samples: Any) -> List[Dict[str, str]]:
    """Accept a TSV path or an inline list of records -> [{sample, r1, r2?, platform?}]."""
    if isinstance(samples, list):
        return samples
    rows: List[Dict[str, str]] = []
    with open(samples) as fh:
        for rec in csv.DictReader(fh, delimiter="\t"):
            rows.append({k: (v or "").strip() for k, v in rec.items()})
    return rows


def _declared_class(platform: str) -> str:
    p = (platform or "illumina").lower()
    if p in _SHORT:
        return "illumina"
    if p in {"ont", "nanopore"}:
        return "ont"
    if p == "pacbio_hifi":
        return "pacbio_hifi"
    if p in {"pacbio_clr", "pacbio"}:
        return "pacbio_clr"
    return p


def _is_long(cls: str) -> bool:
    return cls in {"ont", "pacbio_hifi", "pacbio_clr"}


def _infer_class(median_len: float, est_error: Optional[float]) -> str:
    if median_len <= LONG_MIN_LEN:
        return "illumina"
    if est_error is None:
        return "ont"  # long but no quality (FASTA) -> generic long-noisy
    return "pacbio_hifi" if est_error < HIFI_MAX_ERR else "ont"


# --------------------------------------------------------------------------- #
# Per-file profiling (bounded, local, non-reconstructive)                      #
# --------------------------------------------------------------------------- #
def profile_file(path: str, max_reads: int = 100_000) -> Dict[str, Any]:
    """Aggregate stats over the first ``max_reads`` records. No sequence/ID is retained."""
    fmt = read_format(path)
    lengths: List[int] = []
    gc = bases = 0
    q_err_sum = q_sum = q20 = qbases = 0
    prefix_hashes = set()
    sampled_bytes = 0
    n = 0
    exhausted = True

    with _open(path) as fh:
        if fmt == "fastq":
            while n < max_reads:
                h = fh.readline()
                if not h:
                    break
                s = fh.readline().rstrip("\n")
                fh.readline()  # '+'
                q = fh.readline().rstrip("\n")
                n += 1
                sampled_bytes += len(h) + len(s) + len(q) + 4
                _acc_seq(s, lengths, prefix_hashes)
                u = s.upper()
                gc += u.count("G") + u.count("C")
                bases += len(s)
                for ch in q:
                    Q = ord(ch) - 33
                    q_sum += Q
                    q_err_sum += 10 ** (-Q / 10)
                    if Q >= 20:
                        q20 += 1
                    qbases += 1
            else:
                exhausted = False
        else:  # fasta
            seq_parts: List[str] = []
            started = False
            for line in fh:
                if line.startswith(">"):
                    if started:
                        if _emit_fasta(seq_parts, lengths, prefix_hashes):
                            n += 1
                            s = "".join(seq_parts)
                            sampled_bytes += len(s) + 2
                            gc += s.upper().count("G") + s.upper().count("C")
                            bases += len(s)
                        if n >= max_reads:
                            exhausted = False
                            break
                    seq_parts, started = [], True
                elif started:
                    seq_parts.append(line.rstrip("\n"))
            else:
                if started and n < max_reads and _emit_fasta(seq_parts, lengths, prefix_hashes):
                    n += 1
                    s = "".join(seq_parts)
                    sampled_bytes += len(s) + 2
                    gc += s.upper().count("G") + s.upper().count("C")
                    bases += len(s)

    median_len = _median([float(x) for x in lengths])
    est_error = round(q_err_sum / qbases, 5) if qbases else None
    profile: Dict[str, Any] = {
        "format": fmt,
        "n_sampled": n,
        "read_length": {
            "min": lengths and min(lengths) or 0,
            "median": int(median_len),
            "p90": _percentile(lengths, 0.90),
            "max": lengths and max(lengths) or 0,
        },
        "gc_fraction": round(gc / bases, 3) if bases else None,
        "mean_q": round(q_sum / qbases, 1) if qbases else None,
        "q20_frac": round(q20 / qbases, 3) if qbases else None,
        "est_error": est_error,
        "dup_fraction": round(1 - len(prefix_hashes) / n, 3) if n else None,
        "inferred_platform_class": _infer_class(median_len, est_error),
        "host_fraction": None,  # MVP: reference-dependent, out of scope
        "estimated_bases": _estimate_bases(path, sampled_bytes, n, median_len, bases, exhausted),
    }
    return profile


def _acc_seq(s: str, lengths: List[int], prefix_hashes: set) -> None:
    lengths.append(len(s))
    prefix_hashes.add(hash(s[:PREFIX_BP]))  # store the hash, never the bases


def _emit_fasta(seq_parts: List[str], lengths: List[int], prefix_hashes: set) -> bool:
    s = "".join(seq_parts)
    if not s:
        return False
    _acc_seq(s, lengths, prefix_hashes)
    return True


def _estimate_bases(path: str, sampled_bytes: int, n: int, median_len: float,
                    bases: int, exhausted: bool) -> Optional[float]:
    """Total bases in the file. Exact when fully read; else size-extrapolated (approximate)."""
    if n == 0:
        return 0
    if exhausted:
        return float(bases)
    try:
        size = os.path.getsize(path)
    except OSError:
        return None
    uncomp = size * GZIP_RATIO if is_gzipped(path) else size
    mean_rec_bytes = sampled_bytes / n
    est_reads = uncomp / mean_rec_bytes if mean_rec_bytes else n
    return round(est_reads * median_len)


# --------------------------------------------------------------------------- #
# Reconcile across all samples + build the interview context                   #
# --------------------------------------------------------------------------- #
def reconcile(samples: Dict[str, Dict[str, Any]]) -> Dict[str, Any]:
    inferred = {s: p["inferred_platform_class"] for s, p in samples.items()}
    consensus = next(iter(set(inferred.values()))) if len(set(inferred.values())) == 1 else "mixed"
    medians = [p["read_length"]["median"] for p in samples.values()]
    warnings: List[str] = []
    for s, p in samples.items():
        decl = _declared_class(p.get("declared_platform", "illumina"))
        if _is_long(decl) != _is_long(p["inferred_platform_class"]):
            warnings.append(
                f"sample '{s}': declared {decl} but reads look {p['inferred_platform_class']} "
                "(short/long mismatch) — check the sample sheet."
            )
        if p.get("q20_frac") is not None and p["q20_frac"] < LOW_Q20_FRAC:
            warnings.append(f"sample '{s}': low quality (Q20 fraction {p['q20_frac']}).")
    return {
        "n_samples": len(samples),
        "read_length_median": {
            "min": min(medians) if medians else 0,
            "median": int(_median([float(m) for m in medians])),
            "max": max(medians) if medians else 0,
        },
        "platform_consensus": consensus,
        "warnings": warnings,
    }


def to_context(samples: Dict[str, Dict[str, Any]], project: Dict[str, Any]) -> Dict[str, Any]:
    """Flatten the project view into the promotion-ready context dict (safe reducers)."""
    errs = [p["est_error"] for p in samples.values() if p.get("est_error") is not None]
    q20s = [p["q20_frac"] for p in samples.values() if p.get("q20_frac") is not None]
    bases = [p["estimated_bases"] for p in samples.values() if p.get("estimated_bases")]
    mismatch = any("mismatch" in w for w in project["warnings"])
    return {
        "estimated_bases": max(bases) if bases else None,   # ANY sample deep -> surface asm-coverage
        "platform_class": project["platform_consensus"] if project["platform_consensus"] != "mixed" else None,
        "max_est_error": max(errs) if errs else None,
        "max_host_fraction": 0,                              # MVP: host not measured
        "any_sample_low_q": any(q < LOW_Q20_FRAC for q in q20s),
        "platform_mismatch": mismatch,
        "measured": True,
    }


# --------------------------------------------------------------------------- #
# Entry point                                                                  #
# --------------------------------------------------------------------------- #
def run(samples: Any, max_reads: int = 100_000, max_samples: Optional[int] = None,
        out: Optional[str] = None, assume_yes: bool = False) -> Dict[str, Any]:
    """Probe every sample (subject to consent) and return the report + context dict.

    Consent: ``assume_yes`` records local consent; otherwise the stored choice is used. With
    no consent (or 'off', or non-interactive with nothing stored) this reads nothing and
    returns an advisory stub so the caller falls back to a-priori suggestions.
    """
    decision = consent.set("probe", "local") if assume_yes else consent.get("probe")
    if decision != "local":
        return {"ok": True, "measured": False, "consent": decision,
                "reason": "probe consent not granted; staying advisory "
                          "(run with --yes to allow local read profiling)",
                "context": {"measured": False}}

    rows = load_sheet(samples)
    if max_samples is not None:
        rows = rows[:max_samples]
    profiles: Dict[str, Dict[str, Any]] = {}
    for rec in rows:
        name = rec.get("sample") or rec.get("name")
        r1 = rec.get("r1")
        if not name or not r1 or not os.path.exists(r1):
            continue
        prof = profile_file(r1, max_reads=max_reads)
        prof["declared_platform"] = rec.get("platform", "illumina")
        profiles[name] = prof

    project = reconcile(profiles)
    context = to_context(profiles, project)
    report = {"ok": True, "measured": True, "consent": "local",
              "samples": profiles, "project": project, "context": context}

    if out:
        os.makedirs(out, exist_ok=True)
        import json
        with open(os.path.join(out, "probe.json"), "w") as fh:
            json.dump(report, fh, indent=2)
        with open(os.path.join(out, "probe.md"), "w") as fh:
            fh.write(_render_md(report))
    return report


def _render_md(report: Dict[str, Any]) -> str:
    p = report["project"]
    lines = [f"# Probe report ({p['n_samples']} samples)", "",
             f"- platform consensus: **{p['platform_consensus']}**",
             f"- read-length median across samples: {p['read_length_median']}", ""]
    if p["warnings"]:
        lines += ["## Warnings", ""] + [f"- {w}" for w in p["warnings"]] + [""]
    lines += ["## Per-sample", "", "| sample | platform (declared→inferred) | median len | mean Q | est err | est bases |",
              "|---|---|---|---|---|---|"]
    for s, pr in report["samples"].items():
        lines.append(
            f"| {s} | {pr['declared_platform']}→{pr['inferred_platform_class']} | "
            f"{pr['read_length']['median']} | {pr['mean_q']} | {pr['est_error']} | {pr['estimated_bases']} |"
        )
    return "\n".join(lines) + "\n"
