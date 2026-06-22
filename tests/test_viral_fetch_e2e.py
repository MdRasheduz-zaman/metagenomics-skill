"""Network-gated end-to-end test: classify the noisy ONT fixture reads against the FULL
RefSeq viral database (real NCBI taxonomy), via the prebuilt `fetch-db viral` index.

This is the "sensitivity on noisy long reads + does it find the right viruses among the whole
viral kingdom" check. It is heavier than test_pipeline_e2e (a ~0.6 GB download), so it is
**opt-in** and skips cleanly unless either:
  - METAGX_VIRAL_DB points at an already-fetched RefSeq viral kraken2 index, OR
  - METAGX_E2E_VIRAL_FETCH=1 is set (then the test downloads it with `dbfetch`).
plus kraken2 on PATH. So CI (no network/tools) stays green; run it deliberately:

    export PATH="$HOME/miniconda3/envs/metagx-bio/bin:$PATH"
    METAGX_E2E_VIRAL_FETCH=1 pytest tests/test_viral_fetch_e2e.py -q
    # or, against an index you already have:
    METAGX_VIRAL_DB=/path/to/k2_viral pytest tests/test_viral_fetch_e2e.py -q

Observed truth (30-genome fixture simulated into 1500 ONT reads): ~49% classified at
confidence 0, ~29 distinct species recovered, real ICTV names (Orthoflavivirus denguei /
Bandavirus heartlandense / Alphainfluenzavirus influenzae). The bands below are loose enough
to tolerate RefSeq-viral and kraken2 version drift while still proving real recovery.
"""
from __future__ import annotations

import os
import shutil
from pathlib import Path

import pytest
import yaml

from metagx import dbfetch, runner

REPO = Path(__file__).resolve().parents[1]
FIXTURE_ONT = REPO / "tests" / "fixtures" / "viral" / "ont_reads.fasta"

_HAVE_KRAKEN2 = shutil.which("kraken2") is not None
_VIRAL_DB = os.environ.get("METAGX_VIRAL_DB")
_FETCH = bool(os.environ.get("METAGX_E2E_VIRAL_FETCH"))


def _db_present(path: str | None) -> bool:
    return bool(path) and (Path(path) / "hash.k2d").is_file()


requires_viral = pytest.mark.skipif(
    not (_HAVE_KRAKEN2 and FIXTURE_ONT.is_file() and (_db_present(_VIRAL_DB) or _FETCH)),
    reason="needs kraken2 + the ONT fixture + a RefSeq viral DB (set METAGX_VIRAL_DB to an "
    "index, or METAGX_E2E_VIRAL_FETCH=1 to download ~0.6 GB); skips in CI",
)


@pytest.fixture(scope="session")
def viral_db(tmp_path_factory) -> Path:
    """The full RefSeq viral kraken2 index — a provided one, else fetched once (~0.6 GB)."""
    if _db_present(_VIRAL_DB):
        return Path(_VIRAL_DB)
    db_dir = tmp_path_factory.mktemp("k2_viral")
    res = dbfetch.fetch(name="viral", db_dir=str(db_dir), run=True, force=False)
    assert res.get("ok"), f"viral index fetch failed: {res}"
    assert (db_dir / "hash.k2d").is_file(), f"no hash.k2d after fetch: {list(db_dir.iterdir())}"
    return db_dir


def _parse_kreport(path: Path):
    """(classified_fraction, {species_name: clade_reads})."""
    unclassified = classified = 0
    species = {}
    for line in path.read_text().splitlines():
        c = line.split("\t")
        if len(c) < 6:
            continue
        reads, rank, name = int(c[1]), c[3].strip(), c[5].strip()
        if rank == "U":
            unclassified += reads
        elif rank == "R":
            classified += reads
        elif rank == "S" and reads > 0:
            species[name] = reads
    total = unclassified + classified
    return (classified / total if total else 0.0), species


@requires_viral
def test_ont_reads_recovered_against_full_viral_db(viral_db, tmp_path):
    """Noisy ONT reads classified against the whole RefSeq viral kingdom recover the right
    viruses (real taxonomy), at the ~half-classified rate expected of ONT error rates."""
    from metagx import config_builder as cb

    cfg = cb.build_config(
        project="viralfetch", outdir=str(tmp_path / "out"), threads=4,
        samples=[{"sample": "ont", "r1": str(FIXTURE_ONT), "platform": "ont", "layout": "se"}],
        db={"kraken2": str(viral_db)},
        # abundance off: the prebuilt viral index ships Bracken distributions up to 300mers,
        # not the 1000 an ONT run wants — detection (kraken2) is the sensitivity question here.
        modules={"classify": True, "abundance": False},
        kraken2={"confidence": 0.0})
    cfg_path = tmp_path / "config.yaml"
    cfg_path.write_text(yaml.safe_dump(cfg))

    proc = runner.run(config=str(cfg_path), cores=4)
    assert proc.returncode == 0, (proc.stdout + proc.stderr)[-3000:]

    kreport = next((tmp_path / "out").rglob("*.kreport"))
    frac, species = _parse_kreport(kreport)

    # ~half the noisy ONT reads classify at confidence 0 (the ONT-error signature); loose band.
    assert 0.30 <= frac <= 0.70, f"unexpected classified fraction {frac:.2f}"
    # near-complete recovery of the 30-genome truth set against the full viral kingdom
    assert len(species) >= 18, f"only {len(species)} species recovered: {sorted(species)[:5]}"
    # the right viruses, by real ICTV genus names (robust to species-epithet drift)
    blob = " ".join(species)
    hits = [g for g in ("Orthoflavivirus", "Bandavirus", "Alphainfluenzavirus",
                        "Orthonairovirus", "Lentivirus") if g in blob]
    assert len(hits) >= 2, f"expected known viral genera, recovered: {sorted(species)[:10]}"
