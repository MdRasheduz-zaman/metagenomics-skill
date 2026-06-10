"""Format-aware, seeded read subsampling — no external tools required.

Keeps each record with probability ``fraction`` using a seeded RNG, so runs are
reproducible. Handles FASTA and FASTQ, gzipped or not. Output is written uncompressed in
the input's format. Multi-line FASTA sequences are supported.
"""

from __future__ import annotations

import gzip
import random
from typing import Iterator, Tuple

from .formats import is_gzipped, read_format


def _open(path: str, mode: str = "rt"):
    return gzip.open(path, mode) if is_gzipped(path) else open(path, mode)


def _iter_fasta(fh) -> Iterator[Tuple[str, str]]:
    header, seq = None, []
    for line in fh:
        line = line.rstrip("\n")
        if line.startswith(">"):
            if header is not None:
                yield header, "\n".join(seq)
            header, seq = line, []
        elif header is not None:
            seq.append(line)
    if header is not None:
        yield header, "\n".join(seq)


def _iter_fastq(fh) -> Iterator[Tuple[str, str, str, str]]:
    while True:
        h = fh.readline()
        if not h:
            return
        s = fh.readline()
        plus = fh.readline()
        q = fh.readline()
        yield h.rstrip("\n"), s.rstrip("\n"), plus.rstrip("\n"), q.rstrip("\n")


def subsample(infile: str, outfile: str, fraction: float, seed: int = 42) -> dict:
    """Subsample ``infile`` -> ``outfile``. Returns {kept, total, fraction, seed, format}."""
    if not (0 < fraction <= 1):
        raise ValueError(f"fraction must be in (0, 1], got {fraction}")
    fmt = read_format(infile)
    rng = random.Random(seed)
    total = kept = 0
    with _open(infile) as fin, open(outfile, "wt") as fout:
        if fmt == "fasta":
            for header, seq in _iter_fasta(fin):
                total += 1
                if rng.random() < fraction:
                    kept += 1
                    fout.write(f"{header}\n{seq}\n")
        else:
            for h, s, plus, q in _iter_fastq(fin):
                total += 1
                if rng.random() < fraction:
                    kept += 1
                    fout.write(f"{h}\n{s}\n{plus}\n{q}\n")
    return {"kept": kept, "total": total, "fraction": fraction, "seed": seed, "format": fmt}
