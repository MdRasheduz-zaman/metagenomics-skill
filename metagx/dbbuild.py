"""Build a custom kraken2 + Bracken database from a FASTA of reference genomes.

Self-contained: instead of downloading NCBI taxonomy, it fabricates a minimal taxonomy
(each genome becomes a species directly under root) and tags every sequence with a
``kraken:taxid|`` header, which is all kraken2-build needs. Masking is disabled so no
``dustmasker`` dependency is required.
"""

from __future__ import annotations

import os
import shutil
import subprocess
from typing import Dict, List, Tuple

from .formats import is_gzipped
import gzip

FIRST_TAXID = 1001  # synthetic species taxids start here; root is 1


def _open(path: str):
    return gzip.open(path, "rt") if is_gzipped(path) else open(path, "rt")


def _parse_genomes(genomes) -> List[Tuple[str, str, List[str]]]:
    """Return [(accession, description, [seq_lines])] for each record.

    ``genomes`` may be a single FASTA path or a list of them (a folder of per-genome
    FASTAs), so the synthetic-taxonomy path serves both custom-fasta and custom-folder.
    """
    sources = [genomes] if isinstance(genomes, str) else list(genomes)
    records: List[Tuple[str, str, List[str]]] = []
    for src in sources:
        acc, desc, seq = None, "", []
        with _open(src) as fh:
            for line in fh:
                line = line.rstrip("\n")
                if line.startswith(">"):
                    if acc is not None:
                        records.append((acc, desc, seq))
                    head = line[1:].split(None, 1)
                    acc = head[0]
                    desc = head[1] if len(head) > 1 else acc
                    seq = []
                elif acc is not None:
                    seq.append(line)
        if acc is not None:
            records.append((acc, desc, seq))
    return records


_FASTA_EXT = (".fa", ".fna", ".fasta", ".fa.gz", ".fna.gz", ".fasta.gz")


def _collect_fastas(source: str) -> List[str]:
    """A db.build source -> list of FASTA files (a folder is expanded, a file is wrapped)."""
    if os.path.isdir(source):
        files = sorted(os.path.join(source, f) for f in os.listdir(source)
                       if f.lower().endswith(_FASTA_EXT))
        if not files:
            raise ValueError(f"no FASTA files (*.fa/.fna/.fasta[.gz]) in folder {source}")
        return files
    if os.path.isfile(source):
        return [source]
    raise ValueError(f"db.build source not found: {source}")


def write_library_and_taxonomy(genomes: str, db_dir: str) -> Dict[str, Dict]:
    """Write DB/taxonomy/{nodes,names}.dmp and a tagged library FASTA.

    Returns {accession: {taxid, name}} mapping.
    """
    records = _parse_genomes(genomes)
    if not records:
        raise ValueError(f"No sequences found in {genomes}")

    taxdir = os.path.join(db_dir, "taxonomy")
    os.makedirs(taxdir, exist_ok=True)
    library = os.path.join(db_dir, "custom_library.fasta")

    mapping: Dict[str, Dict] = {}
    nodes = ["1\t|\t1\t|\tno rank\t|\t-\t|"]  # root points to itself
    names = ["1\t|\troot\t|\t\t|\tscientific name\t|"]

    with open(library, "w") as lib:
        for i, (acc, desc, seq) in enumerate(records):
            taxid = FIRST_TAXID + i
            name = desc.replace("\t", " ").strip() or acc
            mapping[acc] = {"taxid": taxid, "name": name}
            nodes.append(f"{taxid}\t|\t1\t|\tspecies\t|\t-\t|")
            names.append(f"{taxid}\t|\t{name}\t|\t\t|\tscientific name\t|")
            # kraken2 reads the taxid from the |kraken:taxid|<n> header token.
            lib.write(f">{acc}|kraken:taxid|{taxid} {name}\n")
            lib.write("\n".join(seq) + "\n")

    with open(os.path.join(taxdir, "nodes.dmp"), "w") as fh:
        fh.write("\n".join(nodes) + "\n")
    with open(os.path.join(taxdir, "names.dmp"), "w") as fh:
        fh.write("\n".join(names) + "\n")
    return mapping


def _have(tool: str) -> bool:
    return shutil.which(tool) is not None


def _usable_cpus() -> int:
    """Number of processors actually available to this process.

    Bracken's ``kmer2read_distr`` hard-fails with "thread count exceeds number of processors"
    when ``-t`` exceeds the online CPU count (kraken2 only warns and reduces). On a 2-core CI
    runner that aborts the whole DB build. Prefer the affinity-aware count (respects cgroup/
    taskset limits) and fall back to ``os.cpu_count()``.
    """
    try:
        return max(1, len(os.sched_getaffinity(0)))  # Linux; honors cgroup/taskset pinning
    except AttributeError:
        return max(1, os.cpu_count() or 1)


def _step_artifacts(name: str, db_dir: str) -> List[str]:
    """The output files a build step must produce to count as successful.

    Used to distinguish a *real* failure from kraken2-build's well-known SIGPIPE quirk:
    on small databases ``build_db`` finishes reading the library and closes the pipe before
    the wrapper's ``cat`` is done, so ``cat`` dies with signal 13, ``xargs`` reports failure,
    and ``kraken2-build`` exits non-zero (64) — even though hash.k2d/opts.k2d/taxo.k2d were
    written correctly. Verifying the artifacts exist is the standard, robust workaround.
    """
    if name == "build":
        return [os.path.join(db_dir, f) for f in ("hash.k2d", "opts.k2d", "taxo.k2d")]
    if name.startswith("bracken-build-"):
        length = name.rsplit("-", 1)[1]
        return [os.path.join(db_dir, f"database{length}mers.kmer_distrib")]
    return []


def _artifacts_present(paths: List[str]) -> bool:
    return bool(paths) and all(os.path.isfile(p) and os.path.getsize(p) > 0 for p in paths)


def build_cat_db(genomes: str, db_dir: str, taxonomy_dir: str, run: bool = True) -> Dict:
    """Build a custom CAT (Contig Annotation Tool) database from reference genomes.

    Predicts proteins (prodigal), maps each to its genome's taxid (same synthetic taxonomy
    used for the kraken2 db), and runs ``CAT_pack prepare``. ``taxonomy_dir`` must contain
    names.dmp + nodes.dmp (e.g. the kraken2 db's taxonomy/ written by build_db).
    """
    os.makedirs(db_dir, exist_ok=True)
    proteins = os.path.join(db_dir, "proteins.faa")
    acc2tax = os.path.join(db_dir, "acc2tax.tsv")
    catdb = os.path.join(db_dir, "catdb")

    # genome accession -> synthetic taxid (FASTA order, matching build_db)
    acc_tax, i = {}, 0
    with _open(genomes) as fh:
        for line in fh:
            if line.startswith(">"):
                i += 1
                acc_tax[line[1:].split()[0]] = FIRST_TAXID - 1 + i

    result = {"db": os.path.abspath(catdb), "n_genomes": len(acc_tax),
              "commands": {"prepare": f"CAT_pack prepare --db_fasta {proteins} "
                           f"--acc2tax {acc2tax} --names {taxonomy_dir}/names.dmp "
                           f"--nodes {taxonomy_dir}/nodes.dmp --db_dir {catdb}"}}

    missing = [t for t in ("prodigal", "CAT_pack") if not _have(t)]
    if not run or missing:
        result["ran"] = False
        if missing:
            result["note"] = f"tools not on PATH: {', '.join(missing)} — not executed"
        return result

    p = subprocess.run(["prodigal", "-i", genomes, "-a", proteins, "-p", "meta", "-q"],
                       capture_output=True, text=True)
    if p.returncode != 0:
        return {**result, "ran": True, "ok": False, "failed_step": "prodigal",
                "tail": (p.stderr or "")[-800:]}
    n = 0
    with open(proteins) as fin, open(acc2tax, "w") as out:
        out.write("accession\taccession.version\ttaxid\tgi\n")
        for line in fin:
            if line.startswith(">"):
                pid = line[1:].split()[0]
                t = acc_tax.get(pid.rsplit("_", 1)[0])
                if t:
                    out.write(f"{pid}\t{pid}\t{t}\t0\n")
                    n += 1
    result["n_proteins"] = n
    p = subprocess.run(["CAT_pack", "prepare", "--db_fasta", proteins, "--acc2tax", acc2tax,
                        "--names", f"{taxonomy_dir}/names.dmp", "--nodes", f"{taxonomy_dir}/nodes.dmp",
                        "--db_dir", catdb], capture_output=True, text=True)
    result["ran"] = True
    result["ok"] = p.returncode == 0
    if p.returncode != 0:
        result["tail"] = ((p.stdout or "") + (p.stderr or ""))[-1000:]
    else:
        result["database_folder"] = os.path.join(catdb, "db")
        result["taxonomy_folder"] = os.path.join(catdb, "tax")
    return result


def build_kaiju_db(genomes: str, db_dir: str, taxonomy_dir: str,
                   threads: int = 4, run: bool = True) -> Dict:
    """Build a custom Kaiju (protein) database from reference genomes — no NCBI download.

    Predicts proteins (prodigal), labels each with its genome's synthetic taxid (Kaiju reads
    the taxid as the protein header), then runs kaiju-mkbwt + kaiju-mkfmi. Copies
    names.dmp/nodes.dmp from ``taxonomy_dir`` (e.g. the kraken2 db's taxonomy/ from build_db)
    so the directory is a drop-in ``db.kaiju`` for the consensus module: it ends up holding
    ``kaiju_db.fmi`` + ``nodes.dmp`` + ``names.dmp``, exactly what rules/consensus.smk expects.
    """
    os.makedirs(db_dir, exist_ok=True)
    raw_proteins = os.path.join(db_dir, "proteins.faa")
    kaiju_proteins = os.path.join(db_dir, "kaiju_proteins.faa")
    prefix = os.path.join(db_dir, "kaiju_db")  # -> kaiju_db.bwt/.sa/.fmi

    # genome accession -> synthetic taxid (FASTA order, matching build_db/build_cat_db)
    acc_tax, i = {}, 0
    with _open(genomes) as fh:
        for line in fh:
            if line.startswith(">"):
                i += 1
                acc_tax[line[1:].split()[0]] = FIRST_TAXID - 1 + i

    result = {"db": os.path.abspath(db_dir), "n_genomes": len(acc_tax),
              "fmi": os.path.abspath(prefix + ".fmi"),
              "commands": {
                  "prodigal": f"prodigal -i {genomes} -a {raw_proteins} -p meta -q",
                  "mkbwt": f"kaiju-mkbwt -n {threads} -a protein -o {prefix} {kaiju_proteins}",
                  "mkfmi": f"kaiju-mkfmi {prefix}"}}

    missing = [t for t in ("prodigal", "kaiju-mkbwt", "kaiju-mkfmi") if not _have(t)]
    if not run or missing:
        result["ran"] = False
        if missing:
            result["note"] = f"tools not on PATH: {', '.join(missing)} — not executed"
        return result

    p = subprocess.run(["prodigal", "-i", genomes, "-a", raw_proteins, "-p", "meta", "-q"],
                       capture_output=True, text=True)
    if p.returncode != 0:
        return {**result, "ran": True, "ok": False, "failed_step": "prodigal",
                "tail": (p.stderr or "")[-800:]}

    # Relabel each protein header to its taxid (Kaiju parses the taxid from the header),
    # and drop prodigal's trailing '*' stop char which kaiju-mkbwt rejects.
    n, total_aa = 0, 0
    with open(raw_proteins) as fin, open(kaiju_proteins, "w") as out:
        seq, taxid = [], None

        def _flush():
            nonlocal total_aa
            if taxid and seq:
                s = "".join(seq).replace("*", "")
                out.write(f">{taxid}\n{s}\n")
                total_aa += len(s)
                return 1
            return 0
        for line in fin:
            if line.startswith(">"):
                n += _flush()
                seq = []
                pid = line[1:].split()[0]               # prodigal: <contig>_<gene>
                taxid = acc_tax.get(pid.rsplit("_", 1)[0])
            else:
                seq.append(line.strip())
        n += _flush()
    result["n_proteins"] = n

    # kaiju-mkbwt estimates its buffer from file size, which underflows for small custom
    # DBs ("Not enough memory allocated"); pass an explicit length in millions (>=1, rounded up).
    import math
    length_mb = max(1, math.ceil(total_aa / 1_000_000))
    result["commands"]["mkbwt"] = (
        f"kaiju-mkbwt -n {threads} -a protein -l {length_mb} -o {prefix} {kaiju_proteins}")
    for step, cmd in (("mkbwt", ["kaiju-mkbwt", "-n", str(threads), "-a", "protein",
                                 "-l", str(length_mb), "-o", prefix, kaiju_proteins]),
                      ("mkfmi", ["kaiju-mkfmi", prefix])):
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0:
            return {**result, "ran": True, "ok": False, "failed_step": step,
                    "tail": ((proc.stdout or "") + (proc.stderr or ""))[-1000:]}

    # copy the taxonomy alongside the index so the dir is a self-contained db.kaiju
    for dmp in ("nodes.dmp", "names.dmp"):
        src = os.path.join(taxonomy_dir, dmp)
        if os.path.isfile(src):
            shutil.copy(src, os.path.join(db_dir, dmp))
    result.update(ran=True, ok=os.path.isfile(prefix + ".fmi"))
    return result


def _execute_steps(steps: List[Tuple[str, List[str]]], db_dir: str) -> Dict:
    """Run kraken2-build/bracken-build steps with the SIGPIPE / thread>core recovery logic.

    Returns a dict to merge into a build result: ``ran``/``ok``[/``failed_step``]/``logs``
    [/``skipped``][/``recovered``][/``note``]. Shared by ``build_db`` (custom synthetic) and
    ``build_database`` (all strategies) so the hard-won recovery rules live in one place.
    """
    logs, skipped, recovered = {}, [], []
    for name, cmd in steps:
        tool = cmd[0]
        if not _have(tool):
            skipped.append(name)
            logs[name] = {"skipped": f"{tool} not on PATH"}
            continue
        proc = subprocess.run(cmd, capture_output=True, text=True)
        logs[name] = {
            "returncode": proc.returncode,
            "tail": ((proc.stdout or "") + (proc.stderr or ""))[-1500:],
        }
        if proc.returncode != 0:
            # A non-zero exit is only a real failure if the step's artifacts are missing.
            # kraken2-build emits a SIGPIPE-driven exit 64 on small DBs while still writing a
            # valid database; trust the artifacts over the exit code.
            artifacts = _step_artifacts(name, db_dir)
            if _artifacts_present(artifacts):
                logs[name]["recovered"] = (
                    f"{tool} exited {proc.returncode}, but {', '.join(os.path.basename(a) for a in artifacts)} "
                    "were produced — treating as success (known SIGPIPE quirk of the build wrapper)."
                )
                recovered.append(name)
                continue
            # Missing artifacts → the kraken2 `--build` genuinely aborted. On low-core CI runners
            # build_db caps OMP threads and its internal `cat | build_db` pipe races so `cat` dies
            # with SIGPIPE *before* step 3 writes any `*.k2d` (seen only on the 2-core GitHub
            # Linux runner; multi-core builds finish cleanly). A single-threaded rebuild removes
            # the race. Retry the build step once with --threads 1 before failing hard.
            if name == "build" and "--threads" in cmd:
                retry_cmd = list(cmd)
                retry_cmd[retry_cmd.index("--threads") + 1] = "1"
                # --threads 1 caps build_db's own threads, but its libgomp regions can still
                # spawn the racing reader; force OMP_NUM_THREADS=1 in the env too so the
                # `cat | build_db` pipe is genuinely serial.
                retry_env = {**os.environ, "OMP_NUM_THREADS": "1"}
                rp = subprocess.run(retry_cmd, capture_output=True, text=True, env=retry_env)
                logs[name]["retry_threads1"] = {
                    "returncode": rp.returncode,
                    "tail": ((rp.stdout or "") + (rp.stderr or ""))[-1500:],
                }
                if rp.returncode == 0 or _artifacts_present(artifacts):
                    logs[name]["recovered"] = (
                        f"{tool} aborted under multithreading (exit {proc.returncode}, no artifacts); "
                        "single-threaded rebuild produced the database — recovered."
                    )
                    recovered.append(name)
                    continue
            return {"ran": True, "ok": False, "failed_step": name, "logs": logs}
    out: Dict = {"ran": True, "ok": True, "logs": logs}
    notes = []
    if recovered:
        out["recovered"] = recovered
        notes.append(
            f"recovered {', '.join(recovered)}: exited non-zero (kraken2-build SIGPIPE quirk on "
            "small DBs — either valid artifacts were written anyway, or a single-threaded rebuild "
            "succeeded); see per-step logs for which."
        )
    if skipped:
        out["skipped"] = skipped
        notes.append(
            f"skipped {', '.join(skipped)} (tool missing); "
            "the kraken2 db is usable but abundance (Bracken) won't run until bracken-build is available"
        )
    if notes:
        out["note"] = " ".join(notes)
    return out


def build_db(
    genomes: str,
    db_dir: str,
    read_length=150,
    threads: int = 4,
    kmer_len: int = 35,
    minimizer_len: int = 31,
    run: bool = True,
) -> Dict:
    """Prepare taxonomy + library, then (optionally) run kraken2-build and bracken-build.

    Returns a dict with the mapping, the planned commands, and per-step results. If the
    tools are missing or run=False, the commands are returned but not executed.
    """
    os.makedirs(db_dir, exist_ok=True)
    mapping = write_library_and_taxonomy(genomes, db_dir)
    library = os.path.join(db_dir, "custom_library.fasta")

    # Never request more threads than there are online CPUs: bracken-build's kmer2read_distr
    # aborts (rc 1, "thread count exceeds number of processors") rather than reducing, which
    # killed the build on the 2-core CI runner. Clamping also keeps kraken2-build from the
    # thread>core mismatch behind its SIGPIPE race.
    requested_threads = threads
    threads = max(1, min(int(threads), _usable_cpus()))

    # Bracken's k-mer distribution is read-length specific; build one per requested length
    # so short- and long-read samples can each use a matching distribution.
    lengths = read_length if isinstance(read_length, (list, tuple)) else [read_length]
    steps = [
        ("add-to-library",
         ["kraken2-build", "--add-to-library", library, "--db", db_dir, "--no-masking"]),
        ("build",
         ["kraken2-build", "--build", "--db", db_dir, "--threads", str(threads),
          "--kmer-len", str(kmer_len), "--minimizer-len", str(minimizer_len)]),
    ]
    for L in lengths:
        steps.append((f"bracken-build-{L}",
                      ["bracken-build", "-d", db_dir, "-t", str(threads),
                       "-k", str(kmer_len), "-l", str(L)]))
    plan = {name: " ".join(cmd) for name, cmd in steps}

    result = {
        "db": os.path.abspath(db_dir),
        "n_sequences": len(mapping),
        "taxids": {acc: m["taxid"] for acc, m in mapping.items()},
        "threads": threads,
        "commands": plan,
    }
    if threads != requested_threads:
        result["note_threads"] = (
            f"requested {requested_threads} threads but only {threads} CPU(s) are available; "
            "clamped (bracken-build/kmer2read_distr aborts when threads exceed online CPUs)."
        )

    # kraken2-build is required to build anything; bracken-build is skippable.
    if not run or not _have("kraken2-build"):
        result["ran"] = False
        if not _have("kraken2-build"):
            result["note"] = "kraken2-build not on PATH — commands not executed"
        return result
    result.update(_execute_steps(steps, db_dir))
    return result


def build_database(
    *,
    db_dir: str,
    strategy: str = "standard",
    taxonomy: str = "real",
    libraries=None,
    source: str = None,
    read_lengths=(150,),
    threads: int = 4,
    kmer_len: int = 35,
    minimizer_len: int = 31,
    minimizer_spaces=None,
    max_db_size=None,
    no_masking=None,
    use_ftp: bool = True,
    run: bool = True,
) -> Dict:
    """Build a kraken2 + Bracken DB per a db.build strategy. Dispatches:

      standard      download NCBI taxonomy + one or more libraries, then build (real taxonomy)
      custom-fasta  one multifasta (synthetic taxonomy, or real if headers carry kraken:taxid|)
      custom-folder a folder of per-genome FASTAs (same taxonomy rules as custom-fasta)
      spike-in      custom genomes ADDED to standard libraries (real taxonomy, required)

    Returns the same shape as ``build_db`` (commands + per-step logs + ran/ok). Masking
    defaults off for synthetic (no dustmasker dep) and on for real builds unless overridden.
    """
    os.makedirs(db_dir, exist_ok=True)
    requested_threads = threads
    threads = max(1, min(int(threads), _usable_cpus()))
    lengths = list(read_lengths) if isinstance(read_lengths, (list, tuple)) else [read_lengths]
    libs = [l.strip() for l in str(libraries or "").split(",") if l.strip()]
    sources = _collect_fastas(source) if source else []
    if no_masking is None:
        no_masking = (taxonomy == "synthetic")
    mask = ["--no-masking"] if no_masking else []
    # NCBI deprecated rsync access to ftp.ncbi.nlm.nih.gov, so kraken2-build's default rsync
    # downloads now fail ("connect refused" on port 873). --use-ftp switches to wget, which
    # works — hence the default. Only the NCBI download steps take it.
    ftp = ["--use-ftp"] if use_ftp else []

    result: Dict = {"db": os.path.abspath(db_dir), "strategy": strategy, "taxonomy": taxonomy,
                    "libraries": libs, "read_lengths": lengths, "threads": threads}
    steps: List[Tuple[str, List[str]]] = []

    # NCBI taxonomy: needed for standard/spike-in and for any real-taxonomy custom build.
    # --skip-maps drops the giant accession2taxid maps; standard libraries carry their own
    # seqid->taxid, and custom sequences are expected to carry kraken:taxid| headers.
    if strategy in {"standard", "spike-in"} or (taxonomy == "real" and sources):
        steps.append(("download-taxonomy",
                      ["kraken2-build", "--download-taxonomy", "--skip-maps", "--db", db_dir] + ftp))
    if strategy in {"standard", "spike-in"}:
        for lib in libs:
            steps.append((f"download-library-{lib}",
                          ["kraken2-build", "--download-library", lib, "--db", db_dir] + mask + ftp))

    if strategy in {"custom-fasta", "custom-folder"} and taxonomy == "synthetic":
        # fabricate a flat taxonomy + a kraken:taxid|-tagged library, then add it
        mapping = write_library_and_taxonomy(sources, db_dir)
        result["n_sequences"] = len(mapping)
        result["taxids"] = {acc: m["taxid"] for acc, m in mapping.items()}
        steps.append(("add-to-library",
                      ["kraken2-build", "--add-to-library",
                       os.path.join(db_dir, "custom_library.fasta"), "--db", db_dir] + mask))
    else:
        # real-taxonomy custom, or spike-in: add each source FASTA as-is (its headers must
        # carry real NCBI taxids, e.g. kraken:taxid|<taxid>, to map into the NCBI taxonomy).
        for i, fa in enumerate(sources):
            steps.append((f"add-to-library-{i}",
                          ["kraken2-build", "--add-to-library", fa, "--db", db_dir] + mask))

    build_cmd = ["kraken2-build", "--build", "--db", db_dir, "--threads", str(threads),
                 "--kmer-len", str(kmer_len), "--minimizer-len", str(minimizer_len)]
    if minimizer_spaces is not None:
        build_cmd += ["--minimizer-spaces", str(minimizer_spaces)]
    if max_db_size:
        build_cmd += ["--max-db-size", str(max_db_size)]
    steps.append(("build", build_cmd))
    for L in lengths:
        steps.append((f"bracken-build-{L}",
                      ["bracken-build", "-d", db_dir, "-t", str(threads),
                       "-k", str(kmer_len), "-l", str(L)]))

    result["commands"] = {name: " ".join(cmd) for name, cmd in steps}
    if threads != requested_threads:
        result["note_threads"] = (
            f"requested {requested_threads} threads but only {threads} CPU(s) are available; "
            "clamped (bracken-build/kmer2read_distr aborts when threads exceed online CPUs)."
        )
    if not run or not _have("kraken2-build"):
        result["ran"] = False
        if not _have("kraken2-build"):
            result["note"] = "kraken2-build not on PATH — commands not executed"
        return result
    result.update(_execute_steps(steps, db_dir))
    return result


def db_is_built(db_dir: str, read_lengths=()) -> bool:
    """True iff a usable kraken2 DB (+ the requested Bracken distributions) already exists,
    so the build step is idempotent and a re-run skips it."""
    core = all(os.path.isfile(os.path.join(db_dir, f)) for f in ("hash.k2d", "opts.k2d", "taxo.k2d"))
    brk = all(os.path.isfile(os.path.join(db_dir, f"database{L}mers.kmer_distrib"))
              for L in read_lengths)
    return core and brk


def write_manifest(db_dir: str, result: Dict) -> str:
    """Record DB provenance next to the index (`.metagx_db.json`) for `metagx report`."""
    import json
    from datetime import datetime, timezone
    path = os.path.join(db_dir, ".metagx_db.json")
    manifest = {k: result.get(k) for k in
                ("strategy", "taxonomy", "libraries", "read_lengths", "threads", "commands")}
    manifest["built_at"] = datetime.now(timezone.utc).isoformat()
    manifest["kraken2_version"] = _probe_version("kraken2")
    with open(path, "w") as fh:
        json.dump(manifest, fh, indent=2)
    return path


def _probe_version(tool: str) -> str:
    try:
        out = subprocess.run([tool, "--version"], capture_output=True, text=True)
        return (out.stdout or out.stderr or "").strip().splitlines()[0] if (out.stdout or out.stderr) else ""
    except (OSError, IndexError):
        return ""
