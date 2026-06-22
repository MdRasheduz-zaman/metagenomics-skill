# Build the kraken2 + Bracken database declared by config['db']['build'], as a pipeline step.
# Idempotent: Snakemake skips the rule once the manifest (.metagx_db.json) exists and is newer
# than any source FASTAs. classify/abundance depend on the manifest via db_ready_input(), so
# the DB is built before — and only once for — the whole run. The read lengths come from
# config (db.build.read_lengths, derived from the sample sheet's platforms by config_builder).


def _db_build_source_inputs(wc):
    """Source FASTAs as inputs so changing them rebuilds; standard builds have no file input."""
    src = (DB_BUILD or {}).get("source")
    if not src:
        return []
    if os.path.isdir(src):
        return [os.path.join(src, f) for f in sorted(os.listdir(src))
                if f.lower().endswith((".fa", ".fna", ".fasta", ".fa.gz", ".fna.gz", ".fasta.gz"))]
    return [src] if os.path.isfile(src) else []


rule build_kraken2_db:
    input:
        _db_build_source_inputs,
    output:
        manifest=DB_MANIFEST,
    threads: THREADS
    params:
        db_dir=DB["kraken2"],
    script:
        "../scripts/build_db.py"
