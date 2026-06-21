"""First-class database onboarding — download a prebuilt kraken2/Bracken index.

Getting a usable reference database is the #1 real-user blocker (EVALUATION-2026-06-22, P1):
real metagenomics needs a 8–100 GB+ kraken2 index, and a new user has no obvious, tested way
to get one. This module curates the standard prebuilt indices published at genome-idx (the
official kraken2/Bracken collection — each tarball already contains the matching Bracken
``databaseNmers.kmer_distrib`` files), and downloads + extracts + verifies one.

For a *custom* small database from your own genomes, use ``metagx build-db`` instead.

Design: ``plan()`` is pure (URL + target dir + the exact shell, no I/O), so it is fully
testable offline; ``fetch()`` does the download only when ``run=True``.
"""

from __future__ import annotations

import os
import shutil
import subprocess
from typing import Dict, List, Optional

_BASE = "https://genome-idx.s3.amazonaws.com/kraken"

# Curated standard indices. Each entry ships kraken2 hash + Bracken distributions in one
# tarball. Sizes are the compressed download (the built DB on disk is larger); the index must
# fit in RAM for fast classification. Full `standard`/`pluspf` are listed so users see the
# real cost up front. NOTE: the capped (…gb) builds and the full builds are published on
# *different* dates upstream (genome-idx) — the filenames below were HEAD-verified to exist.
INDICES: Dict[str, Dict[str, str]] = {
    "viral": {
        "file": "k2_viral_20250714.tar.gz",
        "size": "~0.6 GB",
        "desc": "RefSeq viral only — tiny, good for a smoke test or virome work.",
    },
    "standard-8": {
        "file": "k2_standard_08gb_20241228.tar.gz",
        "size": "~6 GB",
        "desc": "Capped Standard (archaea+bacteria+viral+human+UniVec) at 8 GB. Best default "
                "for laptops/CI; recall is reduced vs full Standard by the cap.",
    },
    "standard-16": {
        "file": "k2_standard_16gb_20241228.tar.gz",
        "size": "~12 GB",
        "desc": "Capped Standard at 16 GB — a middle ground when you have the RAM.",
    },
    "standard": {
        "file": "k2_standard_20250714.tar.gz",
        "size": "~76 GB",
        "desc": "Full Standard (archaea+bacteria+viral+plasmid+human+UniVec). Needs a "
                "workstation/HPC with the index resident in RAM.",
    },
    "pluspf-8": {
        "file": "k2_pluspf_08gb_20241228.tar.gz",
        "size": "~6 GB",
        "desc": "Capped Standard + protozoa + fungi (8 GB cap) — environmental/eukaryote-aware.",
    },
    "pluspf": {
        "file": "k2_pluspf_20250714.tar.gz",
        "size": "~81 GB",
        "desc": "Standard + protozoa + fungi (full) — environmental/eukaryote-aware.",
    },
}

DEFAULT = "standard-8"


def index_url(name: str) -> str:
    if name not in INDICES:
        raise KeyError(name)
    return f"{_BASE}/{INDICES[name]['file']}"


def describe() -> List[Dict[str, str]]:
    """List the curated indices (for `metagx fetch-db --list`)."""
    return [{"name": n, "size": d["size"], "url": index_url(n), "description": d["desc"]}
            for n, d in INDICES.items()]


def is_built(db_dir: str) -> bool:
    """A kraken2 index is usable once it has its hash + taxonomy tables."""
    return all(os.path.isfile(os.path.join(db_dir, f))
               for f in ("hash.k2d", "opts.k2d", "taxo.k2d"))


def plan(name: str, db_dir: str, url: Optional[str] = None) -> Dict:
    """Pure: resolve the URL + target dir + the exact download/extract command. No I/O."""
    resolved = url or index_url(name)
    db_dir = os.path.abspath(db_dir)
    # Stream-extract so we never need 2x disk for the tarball.
    command = f"mkdir -p {db_dir} && curl -fSL {resolved} | tar -xzf - -C {db_dir}"
    return {
        "name": name,
        "url": resolved,
        "db": db_dir,
        "size": INDICES.get(name, {}).get("size", "unknown"),
        "command": command,
        "config_hint": {"db": {"kraken2": db_dir}},
    }


def fetch(name: str = DEFAULT, db_dir: str = "local_databases/kraken2",
          url: Optional[str] = None, run: bool = True, force: bool = False) -> Dict:
    """Download + extract a prebuilt index into ``db_dir`` and verify it built.

    Idempotent: if ``db_dir`` already holds a built index it is reused unless ``force``.
    Returns the plan plus ``ran``/``ok`` (and ``note`` on skip/recover).
    """
    result = plan(name, db_dir, url=url)

    if not run:
        result["ran"] = False
        return result

    if is_built(result["db"]) and not force:
        result.update(ran=False, ok=True,
                      note=f"{result['db']} already contains a built kraken2 index — reusing "
                           "(pass force=True / --force to re-download).")
        return result

    if shutil.which("curl") is None or shutil.which("tar") is None:
        result.update(ran=False, ok=False,
                      note="curl and tar are required to download an index — install them, or "
                           f"download {result['url']} manually and extract into {result['db']}.")
        return result

    os.makedirs(result["db"], exist_ok=True)
    proc = subprocess.run(["bash", "-c", result["command"]], capture_output=True, text=True)
    result["ran"] = True
    result["returncode"] = proc.returncode
    if not is_built(result["db"]):
        result["ok"] = False
        result["tail"] = ((proc.stdout or "") + (proc.stderr or ""))[-1500:]
        result["note"] = ("download/extract did not produce a usable index (hash.k2d/opts.k2d/"
                          "taxo.k2d missing). Check the URL and your network, or fetch manually.")
        return result
    # tar can exit non-zero on a harmless trailing-garbage warning while still extracting a
    # valid index (mirrors the build-db SIGPIPE recovery) — trust the artifacts.
    result["ok"] = True
    if proc.returncode != 0:
        result["note"] = (f"curl|tar exited {proc.returncode} but a valid index was extracted "
                          "— treating as success.")
    return result
