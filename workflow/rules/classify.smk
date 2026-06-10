# kraken2 classification, run once per (sample, sweep value). The command line is
# assembled from the kraken2 registry: base params from config["kraken2"] plus the
# swept parameter for this run, plus workflow-managed io/db/threads flags.

def _kraken_cmd(wc, output, threads):
    reads = reads_for_classify(wc.sample)
    paired = is_paired(wc.sample)
    gz = all(str(r).endswith(".gz") for r in reads)

    base = dict(config.get("kraken2", {}))
    base[SWEEP_PARAM] = LABEL_TO_VALUE[wc.label]
    # FASTA reads have no quality scores; drop the FASTQ-only quality filter.
    if all(formats.read_format(str(r)) == "fasta" for r in reads):
        base.pop("minimum_base_quality", None)

    managed = {
        "db": DB["kraken2"],
        "threads": threads,
        "report": output.report,
        "output": output.kraken,
        "paired": paired,
        "gzip_compressed": gz,
    }
    args = registry.render_args("kraken2", base, managed=managed)
    return "kraken2 " + " ".join(args) + " " + " ".join(str(r) for r in reads)


rule kraken2:
    input:
        reads=lambda wc: reads_for_classify(wc.sample),
    output:
        report=f"{OUT}/kraken2/{{sample}}.{{label}}.kreport",
        # per-read calls — kept (not temp) so reconcile can score read-vs-contig accuracy
        kraken=f"{OUT}/kraken2/{{sample}}.{{label}}.kraken",
    threads: THREADS
    params:
        cmd=lambda wc, output, threads: _kraken_cmd(wc, output, threads),
    shell:
        "{params.cmd}"


# Per-sample comparison matrix across the sweep (the "k-dense" view).
rule classify_matrix:
    input:
        reports=lambda wc: expand(
            f"{OUT}/kraken2/{wc.sample}.{{label}}.kreport", label=SWEEP_LABELS
        ),
    output:
        json=f"{OUT}/summary/{{sample}}.matrix.json",
        png=f"{OUT}/summary/{{sample}}.heatmap.png",
    params:
        labels=SWEEP_LABELS,
        sweep_param=SWEEP_PARAM,
        sample=lambda wc: wc.sample,
    script:
        "../scripts/parse_matrix.py"
