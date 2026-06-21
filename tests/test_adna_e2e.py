"""End-to-end test of the ancient-DNA (aDNA) authentication path.

The damage module's science is: map reads to a reference, quantify post-mortem cytosine
deamination with mapDamage2 (elevated C→T at 5' read ends, G→A at 3' ends), and emit an
authenticity verdict (`workflow/scripts/damage_authenticate.py`). The full Snakemake module
maps to a *de novo assembly*, but MEGAHIT segfaults on this Apple-Silicon/Rosetta env, so here
we exercise the aDNA-specific steps directly against the viral reference genomes — which is
exactly what the workflow's mapdamage→authenticate rules do, minus the (shared, arch-broken)
assembly step.

We simulate two libraries from the 30 viral genomes:
  * a damaged library (terminal C→T / G→A deamination applied) -> must authenticate as aDNA,
  * an undamaged control                                       -> must read as modern.

Skips cleanly when mapDamage/minimap2/samtools are absent (e.g. CI).
"""
from __future__ import annotations

import importlib.util
import random
import shutil
import subprocess
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
# Prefer the gitignored full data/ set when present, else fall back to the committed
# 30-genome viral fixture — so this aDNA test actually RUNS in CI (mapDamage is in
# environment.yml) instead of silently skipping for want of gitignored state.
_DATA_GENOMES = REPO / "data" / "genomes.fasta"
_FIXTURE_GENOMES = REPO / "tests" / "fixtures" / "viral" / "genomes.fasta"
GENOMES = _DATA_GENOMES if _DATA_GENOMES.is_file() else _FIXTURE_GENOMES

_HAVE = bool(shutil.which("mapDamage") and shutil.which("minimap2") and shutil.which("samtools"))
requires_mapdamage = pytest.mark.skipif(
    not (_HAVE and GENOMES.is_file()),
    reason="mapDamage/minimap2/samtools not on PATH (skips without the bio stack); "
    "the reference comes from data/ or the committed viral fixture",
)


def _load_authenticate():
    spec = importlib.util.spec_from_file_location(
        "damage_authenticate",
        REPO / "workflow" / "scripts" / "damage_authenticate.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _read_fasta(path: Path):
    name, seq, out = None, [], []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if name:
                out.append((name, "".join(seq)))
            name, seq = line[1:].split()[0], []
        else:
            seq.append(line.strip().upper())
    if name:
        out.append((name, "".join(seq)))
    return out


def _simulate(genomes, dst: Path, n_reads: int, read_len: int, damaged: bool, seed: int):
    """Fragment genomes into short reads, optionally applying terminal deamination.

    Damage model (per read): 5' positions C->T with probability 0.45*exp(-i/3); 3' positions
    G->A with the same decay from the 3' end. This reproduces the double-stranded aDNA
    signature mapDamage looks for. Reads are taken from the forward strand only so the mapped
    orientation matches the deaminated ends.
    """
    rng = random.Random(seed)
    usable = [(nm, s) for nm, s in genomes if len(s) > read_len + 1]
    with open(dst, "w") as fh:
        for i in range(n_reads):
            nm, s = rng.choice(usable)
            start = rng.randint(0, len(s) - read_len)
            read = list(s[start:start + read_len])
            if damaged:
                for j in range(read_len):
                    if read[j] == "C" and rng.random() < 0.45 * (2.718 ** (-j / 3)):
                        read[j] = "T"
                    k = read_len - 1 - j
                    if read[k] == "G" and rng.random() < 0.45 * (2.718 ** (-j / 3)):
                        read[k] = "A"
            seq = "".join(read)
            fh.write(f"@r{i}_{nm}_{start}\n{seq}\n+\n{'I' * read_len}\n")
    return dst


def _map(reads: Path, ref: Path, bam: Path):
    p1 = subprocess.run(["minimap2", "-ax", "sr", str(ref), str(reads)],
                        capture_output=True)
    assert p1.returncode == 0, p1.stderr.decode()[-2000:]
    sort = subprocess.run(["samtools", "sort", "-o", str(bam), "-"], input=p1.stdout,
                          capture_output=True)
    assert sort.returncode == 0, sort.stderr.decode()[-2000:]
    subprocess.run(["samtools", "index", str(bam)], check=True)


def _mapdamage(bam: Path, ref: Path, outdir: Path):
    # --no-stats: just the terminal-frequency tables (skips the R/rpy2 Bayesian model).
    p = subprocess.run(["mapDamage", "-i", str(bam), "-r", str(ref),
                        "-d", str(outdir), "--no-stats"], capture_output=True)
    assert p.returncode == 0, p.stderr.decode()[-2500:]


@requires_mapdamage
@pytest.mark.parametrize("damaged,expect", [(True, True), (False, False)])
def test_adna_damage_authentication(tmp_path, damaged, expect):
    da = _load_authenticate()
    genomes = _read_fasta(GENOMES)
    ref = GENOMES

    reads = _simulate(genomes, tmp_path / "reads.fastq", n_reads=4000, read_len=50,
                      damaged=damaged, seed=7)
    bam = tmp_path / "aln.sorted.bam"
    _map(reads, ref, bam)

    outdir = tmp_path / "mapdamage"
    _mapdamage(bam, ref, outdir)
    ct5 = outdir / "5pCtoT_freq.txt"
    ga3 = outdir / "3pGtoA_freq.txt"
    assert ct5.is_file() and ga3.is_file(), "mapDamage did not emit terminal frequency tables"

    result = da.run(str(ct5), str(ga3), sample="adna_test", threshold=0.05,
                    out_json=str(tmp_path / "auth.json"))

    assert result["damage_present"] is expect, (
        f"damaged={damaged}: got {result}"
    )
    if damaged:
        # the terminal signal should be clearly elevated, not marginal
        assert result["ct_5prime_pos1"] > 0.10
        assert result["ga_3prime_pos1"] > 0.10
