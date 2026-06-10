# Second-classifier consensus (Tier 2), gated by modules.classify_consensus (requires
# modules.classify). Runs an independent classifier — MetaPhlAn (markers) or Kaiju (protein) —
# alongside kraken2 and cross-checks them at the species level. Agreement = high confidence;
# kraken2-only taxa flag DB-completeness false positives. WGS-only. Run with --use-conda.

CONSENSUS_CLASSIFIER = str(config.get("consensus", {}).get("classifier", "metaphlan")).lower()


def second_profile(sample):
    if CONSENSUS_CLASSIFIER == "kaiju":
        return f"{OUT}/consensus/{sample}.kaiju.species.tsv"
    return f"{OUT}/consensus/{sample}.metaphlan.tsv"


# ----------------------- marker-gene classifier: MetaPhlAn -----------------------
rule metaphlan:
    wildcard_constraints:
        sample=_alt(WGS_SAMPLES),
    input:
        reads=lambda wc: reads_for_classify(wc.sample),
    output:
        tsv=f"{OUT}/consensus/{{sample}}.metaphlan.tsv",
    threads: THREADS
    conda:
        "../envs/functional.yaml"
    params:
        reads=lambda wc: ",".join(str(r) for r in reads_for_classify(wc.sample)),
        args=lambda wc, output, threads: " ".join(registry.render_args(
            "metaphlan", config.get("metaphlan", {}),
            managed={"input_type": "fastq", "nproc": threads, "output_file": output.tsv,
                     "bowtie2out": f"{OUT}/consensus/{wc.sample}.bowtie2.bz2",
                     "bowtie2db": DB.get("metaphlan", "")})),
    shell:
        "metaphlan {params.reads} {params.args}"


# ----------------------- protein classifier: Kaiju -----------------------
def _kaiju_cmd(wc, output, threads):
    reads = reads_for_classify(wc.sample)
    kdb = DB.get("kaiju", "")
    nodes = f"{kdb}/nodes.dmp" if kdb else "nodes.dmp"
    names = f"{kdb}/names.dmp" if kdb else "names.dmp"
    raw = f"{OUT}/consensus/{wc.sample}.kaiju.out"
    managed = {"nodes": nodes, "fmi": f"{kdb}/kaiju_db.fmi" if kdb else "kaiju_db.fmi",
               "input1": reads[0], "output": raw, "threads": threads}
    if is_paired(wc.sample):
        managed["input2"] = reads[1]
    args = registry.render_args("kaiju", config.get("kaiju", {}), managed=managed)
    return ("kaiju " + " ".join(args) +
            f" && kaiju2table -t {nodes} -n {names} -r species -o {output.tsv} {raw}")


rule kaiju:
    wildcard_constraints:
        sample=_alt(WGS_SAMPLES),
    input:
        reads=lambda wc: reads_for_classify(wc.sample),
    output:
        tsv=f"{OUT}/consensus/{{sample}}.kaiju.species.tsv",
    threads: THREADS
    conda:
        "../envs/kaiju.yaml"
    params:
        cmd=lambda wc, output, threads: _kaiju_cmd(wc, output, threads),
    shell:
        "{params.cmd}"


# ----------------------- agreement: kraken2 vs the chosen second classifier -----------------------
rule classifier_consensus:
    wildcard_constraints:
        sample=_alt(WGS_SAMPLES),
    input:
        kraken=f"{OUT}/kraken2/{{sample}}.{READ_LABEL}.kreport",
        second=lambda wc: second_profile(wc.sample),
    output:
        json=f"{OUT}/consensus/{{sample}}.consensus.json",
    params:
        classifier=CONSENSUS_CLASSIFIER,
        sample=lambda wc: wc.sample,
    script:
        "../scripts/classifier_consensus.py"
