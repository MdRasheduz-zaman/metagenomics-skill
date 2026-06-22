"""Per-tool reference-database provisioning for the domain/functional modules.

Unlike kraken2/Bracken (one coherent `db.build` — see dbbuild.py), each of these tools ships
its *own* downloader with its own layout, so this is a thin, per-tool fetch layer rather than a
unified builder. Each entry knows: the canonical download command, a presence marker (so the
fetch is idempotent — skipped if the DB is already there), the size, and which module/flag
makes the DB necessary (so `metagx doctor` only demands it when the run will actually use it).

Dependency-light (stdlib only), so it imports cleanly in the CLI, doctor, and the workflow.
"""
from __future__ import annotations

import glob
import os
import shutil
import subprocess
from typing import Callable, Dict, List, Optional

# Each spec: how to fetch the DB, how to tell it's already there, and how big it is.
#   cmd      : db_dir -> argv (the tool's canonical downloader)
#   tool     : the binary that must be on PATH to fetch
#   markers  : globs (relative to db_dir) that exist once the DB is provisioned
#   env      : extra environment for the download command (e.g. GTDB-Tk's data path)
#   size     : rough on-disk size, surfaced to the user (they are often space-constrained)
#   manual   : True => no clean CLI downloader; we can only point the user at the docs
DBSpec = Dict[str, object]

SPECS: Dict[str, DBSpec] = {
    "genomad": {
        "tool": "genomad",
        "cmd": lambda d: ["genomad", "download-database", d],
        "markers": ["genomad_db", "genomad_db/genomad_db.source"],
        "size": "~1.5 GB",
        "needed_by": "domain_taxonomy + viral domain",
    },
    "checkv": {
        "tool": "checkv",
        "cmd": lambda d: ["checkv", "download_database", d],
        "markers": ["checkv-db-*", "checkv-db-*/genome_db/*.dmnd"],
        "size": "~6 GB (extracted; ~0.5 GB download)",  # verified: checkv-db-v1.5 -> 6.4 GB
        "needed_by": "domain_taxonomy + viral domain",
    },
    "checkm2": {
        "tool": "checkm2",
        "cmd": lambda d: ["checkm2", "database", "--download", "--path", d],
        "markers": ["*.dmnd", "**/*.dmnd", "CheckM2_database/*.dmnd"],
        "size": "~3 GB",
        "needed_by": "domain_taxonomy + prokaryote domain",
    },
    "gtdbtk": {
        "tool": "download-db.sh",          # the script GTDB-Tk ships; reads GTDBTK_DATA_PATH
        "cmd": lambda d: ["download-db.sh"],
        "env": lambda d: {"GTDBTK_DATA_PATH": d},
        "markers": ["taxonomy", "markers", "*/taxonomy", "release*"],
        "size": "~110 GB (large — check disk first)",
        "needed_by": "domain_taxonomy + prokaryote domain",
    },
    "bakta": {
        "tool": "bakta_db",
        # default to the LIGHT db (~1.5 GB) over full (~30 GB) — most users + the space-tight case
        "cmd": lambda d: ["bakta_db", "download", "--output", d, "--type", "light"],
        "markers": ["db-light/version.json", "db*/version.json", "**/*.dmnd"],
        "size": "~1.5 GB (light) / ~30 GB (full)",
        "needed_by": "functional.annotation",
    },
    "amrfinderplus": {
        "tool": "amrfinder_update",
        "cmd": lambda d: ["amrfinder_update", "--force_update", "--database", d],
        "markers": ["**/AMR.LIB", "**/AMRProt", "latest/AMR.LIB"],
        "size": "~0.2 GB",
        "needed_by": "functional.amr",
        # the amrfinder rule only runs when db.amrfinderplus is set, so a missing DB is a
        # graceful skip (abricate still covers AMR), not a fatal mid-run crash.
        "self_gates": True,
    },
    "antismash": {
        "tool": "download-antismash-databases",
        "cmd": lambda d: ["download-antismash-databases", "--database-dir", d],
        "markers": ["pfam", "clusterblast", "*/pfam*"],
        "size": "~9 GB",
        "needed_by": "bgc",
    },
    "humann_nucleotide": {
        "tool": "humann_databases",
        "cmd": lambda d: ["humann_databases", "--download", "chocophlan", "full", d, "--update-config", "no"],
        "markers": ["chocophlan*", "**/*.v*.bz2", "**/chocophlan*"],
        "size": "~16 GB",
        "needed_by": "functional.pathways",
    },
    "humann_protein": {
        "tool": "humann_databases",
        "cmd": lambda d: ["humann_databases", "--download", "uniref", "uniref90_diamond", d, "--update-config", "no"],
        "markers": ["uniref*", "**/uniref90*.dmnd"],
        "size": "~20 GB",
        "needed_by": "functional.pathways",
    },
    "eggnog": {
        "tool": "download_eggnog_data.py",
        "cmd": lambda d: ["download_eggnog_data.py", "-y", "--data_dir", d],
        "markers": ["eggnog.db", "**/eggnog.db"],
        "size": "~50 GB",
        "needed_by": "functional.annotation",
    },
    "metaphlan": {
        "tool": "metaphlan",
        "cmd": lambda d: ["metaphlan", "--install", "--bowtie2db", d],
        "markers": ["mpa_*.pkl", "*.pkl", "**/mpa_*"],
        "size": "~25 GB",
        "needed_by": "classify_consensus (when MetaPhlAn is the chosen consensus tool)",
    },
    # --- no clean CLI downloader: gate-only (doctor reports + points at the docs) ----------
    "eukcc": {
        "tool": None,
        "manual": True,
        "docs": "http://eukcc.readthedocs.io (download the EukCC2 DB tarball, set db.eukcc)",
        "markers": ["**/*.dmnd", "**/refpkg", "eukcc2_db*"],
        "size": "~12 GB",
        "needed_by": "domain_taxonomy + eukaryote domain",
    },
    "emu": {
        "tool": None,
        "manual": True,
        "docs": "https://github.com/treangenlab/emu#emu-database (download the default 16S DB, set db.emu)",
        "markers": ["species_taxid.fasta", "**/species_taxid.fasta", "taxonomy.tsv"],
        "size": "~0.7 GB",
        "needed_by": "amplicon (long-read 16S)",
    },
    "blast": {
        # NCBI nt is ~200 GB (like GTDB-Tk: never auto-fetch on a dev box). No clean
        # per-dir CLI downloader (update_blastdb.pl writes to cwd), so this is gate-only;
        # for validation, a small custom DB from makeblastdb or `validate.remote` is the
        # laptop-friendly path. db.blast may be a directory OR a DB-name prefix (e.g. .../nt).
        "tool": None,
        "manual": True,
        "docs": ("get a BLAST+ nucleotide DB: `update_blastdb.pl --decompress nt` (~200 GB, "
                 "needs BLASTDB set) or build a custom one with `makeblastdb -dbtype nucl -in "
                 "refs.fasta -out <dir>/mydb`; or set validate.remote: true to search NCBI"),
        "markers": ["*.nin", "*.nal", "*.ndb", "**/*.nin", "**/*.nal"],
        "prefix": True,
        "prefix_suffixes": ["*.nin", "*.nal", "*.ndb"],
        "size": "~200 GB (nt) / tiny (custom makeblastdb)",
        "needed_by": "validate (BLAST cross-check of classifier calls)",
    },
}


def _have(tool: str) -> bool:
    return shutil.which(tool) is not None


def is_provisioned(tool: str, db_dir: Optional[str]) -> bool:
    """True if the tool's DB already exists at db_dir (any presence marker matches).

    Most tools store their DB in a directory. Some (BLAST) take a DB-name *prefix* (e.g.
    ``.../nt`` for files ``nt.nin``, ``nt.nal``); for those, ``spec['prefix']`` enables a
    sibling-glob check so a valid prefix path isn't reported as missing.
    """
    if not db_dir:
        return False
    spec = SPECS.get(tool)
    if os.path.isdir(db_dir):
        if not spec:
            return bool(os.listdir(db_dir))  # unknown tool: any non-empty dir counts
        for pat in spec["markers"]:           # type: ignore[index]
            if glob.glob(os.path.join(db_dir, pat), recursive=True):
                return True
    # prefix case: db_dir names a DB, not a directory (e.g. BLAST's /path/nt)
    if spec and spec.get("prefix"):
        for suf in spec.get("prefix_suffixes", []):   # type: ignore[union-attr]
            if glob.glob(db_dir + suf):
                return True
    return False


def provision(tool: str, db_dir: str, run: bool = True, force: bool = False) -> Dict:
    """Fetch a tool's DB into db_dir via its canonical downloader. Idempotent: returns early
    (skipped) when the DB is already present unless force=True. Returns a result dict with the
    planned command and outcome (mirrors dbbuild's shape)."""
    spec = SPECS.get(tool)
    if not spec:
        return {"tool": tool, "ok": False, "error": f"no provisioner for '{tool}'; "
                f"known: {sorted(SPECS)}"}
    os.makedirs(db_dir, exist_ok=True)
    if spec.get("manual"):  # no clean CLI downloader — point at the docs (still idempotent-check)
        present = is_provisioned(tool, db_dir)
        return {"tool": tool, "db": os.path.abspath(db_dir), "size": spec.get("size"),
                "ran": False, "ok": present, "manual": True,
                "note": ("already present" if present
                         else f"no automatic downloader for {tool}; get it manually: {spec.get('docs')}")}
    cmd: List[str] = spec["cmd"](db_dir)                       # type: ignore[operator]
    result: Dict = {"tool": tool, "db": os.path.abspath(db_dir), "size": spec.get("size"),
                    "command": " ".join(cmd)}

    if not force and is_provisioned(tool, db_dir):
        result.update(ran=False, ok=True, skipped="already present")
        return result
    binary = spec["tool"]                                      # type: ignore[index]
    if not run or not _have(binary):
        result["ran"] = False
        if not _have(binary):
            result["note"] = f"{binary} not on PATH — command not executed (use --use-conda)"
        return result

    env = dict(os.environ)
    if spec.get("env"):
        env.update(spec["env"](db_dir))                       # type: ignore[operator]
    proc = subprocess.run(cmd, capture_output=True, text=True, env=env)
    result.update(ran=True, ok=(proc.returncode == 0 and is_provisioned(tool, db_dir)),
                  returncode=proc.returncode,
                  tail=((proc.stdout or "") + (proc.stderr or ""))[-1500:])
    return result


def fetch_command(tool: str, db_dir: str = "<dir>") -> str:
    """The canonical downloader command string, for doctor remedies and SKILL docs."""
    spec = SPECS.get(tool)
    if not spec:
        return f"(no metagx provisioner for {tool})"
    if spec.get("manual"):
        return f"(manual download — {spec.get('docs')})"
    pre = ""
    if spec.get("env"):
        pre = " ".join(f"{k}={v}" for k, v in spec["env"](db_dir).items()) + " "  # type: ignore[operator]
    return pre + " ".join(spec["cmd"](db_dir))                # type: ignore[operator]


def needed_dbs(cfg: Dict) -> Dict[str, str]:
    """Map {tool: db_path_key} for every module DB this config's enabled modules will use.
    Drives the doctor presence-gate: a needed DB that's absent fails fast with the fix."""
    mods = cfg.get("modules", {})
    domains = [str(d).lower() for d in cfg.get("domains", [])]
    func = cfg.get("functional", {})
    need: Dict[str, str] = {}
    if mods.get("domain_taxonomy"):
        if "viral" in domains:
            need.update(genomad="genomad", checkv="checkv")
        if "prokaryote" in domains:
            need.update(checkm2="checkm2", gtdbtk="gtdbtk")
        if "eukaryote" in domains:
            need["eukcc"] = "eukcc"
    if mods.get("functional"):
        if func.get("annotation"):
            need.update(bakta="bakta", eggnog="eggnog")
        if func.get("amr"):
            need["amrfinderplus"] = "amrfinderplus"   # note: the rule also self-gates on the db
        if func.get("pathways"):
            need.update(humann_nucleotide="humann_nucleotide", humann_protein="humann_protein")
    if mods.get("bgc"):
        need["antismash"] = "antismash"
    # validate BLASTs classifier calls — needs a local BLAST DB unless searching NCBI remotely.
    if mods.get("validate") and not (cfg.get("validate", {}) or {}).get("remote"):
        need["blast"] = "blast"
    # emu (long-read amplicon only) and metaphlan (an *optional* consensus tool) have
    # platform/sub-tool nuances, so they're fetchable (SPECS) but not auto-required here to
    # avoid false-positive doctor failures — leave them to the explicit path / db.provision.
    return need
