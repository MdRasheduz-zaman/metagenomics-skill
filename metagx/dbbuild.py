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


def _parse_genomes(genomes: str) -> List[Tuple[str, str, List[str]]]:
    """Return [(accession, description, [seq_lines])] for each record."""
    records: List[Tuple[str, str, List[str]]] = []
    acc, desc, seq = None, "", []
    with _open(genomes) as fh:
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
        "commands": plan,
    }

    # kraken2-build is required to build anything; bracken-build is skippable.
    if not run or not _have("kraken2-build"):
        result["ran"] = False
        if not _have("kraken2-build"):
            result["note"] = "kraken2-build not on PATH — commands not executed"
        return result

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
            result.update(ran=True, ok=False, failed_step=name, logs=logs)
            return result
    result.update(ran=True, ok=True, logs=logs)
    notes = []
    if recovered:
        result["recovered"] = recovered
        notes.append(
            f"{', '.join(recovered)} exited non-zero but produced valid artifacts "
            "(kraken2-build SIGPIPE quirk on small DBs) — recovered."
        )
    if skipped:
        result["skipped"] = skipped
        notes.append(
            f"skipped {', '.join(skipped)} (tool missing); "
            "the kraken2 db is usable but abundance (Bracken) won't run until bracken-build is available"
        )
    if notes:
        result["note"] = " ".join(notes)
    return result
