"""Sequence file format helpers, shared by the workflow, CLI, and config validation.

Lets the pipeline handle FASTA (no quality scores) as gracefully as FASTQ: detection is
extension-first with a content sniff fallback, so a `.fasta` of reads classifies fine and
QC (a FASTQ-only step) is skipped automatically for it.
"""

from __future__ import annotations

import csv
import gzip
import os
from typing import Dict, List

FASTQ_EXTS = (".fastq", ".fq")
FASTA_EXTS = (".fasta", ".fa", ".fna", ".fas")


def read_tsv_dicts(path: str) -> List[Dict[str, str]]:
    """Read a TSV into a list of ``{column: value}`` rows, tolerant of real-world sheets.

    The target audience edits sample sheets in Excel, which commonly saves UTF-8 **with a BOM**;
    a naive ``open`` then makes the first header ``"\\ufeffsample"``, so ``row["sample"]`` raises a
    cryptic ``KeyError``. Opening with ``utf-8-sig`` strips the BOM, and we also normalize
    surrounding whitespace on the header names (``" sample "`` -> ``"sample"``). Values are left
    as-is; callers strip where they need to. Every sample-sheet reader routes through here so the
    tolerance can't drift between the CLI, the config builder, the probe, and the workflow.
    """
    with open(path, encoding="utf-8-sig", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        return [{(k.strip() if isinstance(k, str) else k): v for k, v in row.items()}
                for row in reader]


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


def kreport_row(line: str):
    """Parse one kraken2 report line, robust to ``--report-minimizer-data``.

    A standard kraken2 report has 6 tab columns:
        pct, clade_reads, taxon_reads, rank, taxid, name
    ``--report-minimizer-data`` inserts two columns (distinct minimizers, distinct k-mers)
    *between* ``taxon_reads`` and ``rank``, giving 8 columns. The leading three and the
    trailing three (rank, taxid, name) are stable; only the middle grows — so we index the
    rank/taxid/name from the END. Every kreport consumer should use this rather than fixed
    indices, or it silently misreads the rank/name when minimizer reporting is on.

    Returns a dict (``pct, clade_reads, taxon_reads, rank, taxid, name``) or ``None`` for
    malformed / too-short lines. ``name`` keeps its leading indentation (2 spaces per rank)
    because Krona derives lineage depth from it.
    """
    c = line.rstrip("\n").split("\t")
    if len(c) < 6:
        return None
    return {
        "pct": c[0],
        "clade_reads": c[1],
        "taxon_reads": c[2],
        "rank": c[-3].strip(),
        "taxid": c[-2].strip(),
        "name": c[-1],
    }
