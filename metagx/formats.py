"""Sequence file format helpers, shared by the workflow, CLI, and config validation.

Lets the pipeline handle FASTA (no quality scores) as gracefully as FASTQ: detection is
extension-first with a content sniff fallback, so a `.fasta` of reads classifies fine and
QC (a FASTQ-only step) is skipped automatically for it.
"""

from __future__ import annotations

import gzip
import os

FASTQ_EXTS = (".fastq", ".fq")
FASTA_EXTS = (".fasta", ".fa", ".fna", ".fas")


def is_gzipped(path: str) -> bool:
    return str(path).endswith(".gz")


def _strip_gz(path: str) -> str:
    return path[:-3] if path.endswith(".gz") else path


def read_format(path: str) -> str:
    """Return 'fastq' or 'fasta'. Extension-first, with a first-byte sniff fallback."""
    base = _strip_gz(str(path)).lower()
    if base.endswith(FASTQ_EXTS):
        return "fastq"
    if base.endswith(FASTA_EXTS):
        return "fasta"
    # Unknown extension: sniff the first non-empty character.
    try:
        opener = gzip.open if is_gzipped(path) else open
        with opener(path, "rt") as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                if line.startswith("@"):
                    return "fastq"
                if line.startswith(">"):
                    return "fasta"
                break
    except OSError:
        pass
    return "fastq"  # safe default


def is_fasta(path: str) -> bool:
    return read_format(path) == "fasta"


def canonical_ext(path: str) -> str:
    """Canonical extension to use for a derived (e.g. subsampled) file of the same format."""
    return ".fasta" if read_format(path) == "fasta" else ".fastq"


def estimate_read_length(path: str, sample: int = 5000) -> dict:
    """Estimate read-length stats from the first ``sample`` records (FASTA or FASTQ, ±gz).

    Used to pick the Bracken build/run length: Bracken's k-mer distribution is read-length
    specific, so the value must track the data — typical ~150 for Illumina, but the *median*
    for long reads (often 1–10 kb). Returns {median, mean, min, max, n}.
    """
    lengths = []
    opener = gzip.open if is_gzipped(path) else open
    fmt = read_format(path)
    with opener(path, "rt") as fh:
        if fmt == "fastq":
            i = 0
            for line in fh:
                i += 1
                if i % 4 == 2:
                    lengths.append(len(line.strip()))
                    if len(lengths) >= sample:
                        break
        else:
            seq, started = 0, False
            for line in fh:
                if line.startswith(">"):
                    if started:
                        lengths.append(seq)
                        if len(lengths) >= sample:
                            break
                    started, seq = True, 0
                else:
                    seq += len(line.strip())
            if started and len(lengths) < sample and seq:
                lengths.append(seq)
    if not lengths:
        return {"median": 0, "mean": 0, "min": 0, "max": 0, "n": 0}
    lengths.sort()
    n = len(lengths)
    median = lengths[n // 2] if n % 2 else (lengths[n // 2 - 1] + lengths[n // 2]) // 2
    return {"median": median, "mean": round(sum(lengths) / n, 1),
            "min": lengths[0], "max": lengths[-1], "n": n}
