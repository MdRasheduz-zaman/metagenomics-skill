# Contig classification + reconciliation against read-level results.
#   classify_contigs : kraken2 on the assembled contigs (one call per contig)
#   reconcile        : join contig taxonomy + per-contig coverage (depth, reused from the
#                      mapping rules) + read-level kraken2 calls into a per-taxon table,
#                      a concordance view, and discordance/chimera flags.
# Requires modules.assembly (contigs) and modules.classify (read calls + db). Coverage
# comes from the shared map_to_contigs -> contig_depth rules (no need to enable binning).

def _classify_contigs_cmd(wc, output, threads):
    base = dict(config.get("kraken2", {}))
    base.pop(SWEEP_PARAM, None)      # contigs are not swept; one confident pass
    base["use_names"] = True         # emit "name (taxid)" so reconcile can read both
    managed = {"db": DB["kraken2"], "threads": threads,
               "report": output.report, "output": output.calls}
    args = registry.render_args("kraken2", base, managed=managed)
    return "kraken2 " + " ".join(args) + f" {input_contigs(wc)}"


def input_contigs(wc):
    return f"{OUT}/assembly/{wc.sample}/final.contigs.fa"


# Optional CAT (per-ORF contig taxonomy) cross-check, enabled when db.cat is configured.
CAT_DB = DB.get("cat")


def _cat_named(wc):
    return f"{OUT}/cat/{wc.sample}.named.txt" if CAT_DB else []


rule cat_classify:
    input:
        contigs=lambda wc: input_contigs(wc),
    output:
        named=f"{OUT}/cat/{{sample}}.named.txt",
    threads: THREADS
    params:
        prefix=lambda wc: f"{OUT}/cat/{wc.sample}",
        c2c=lambda wc: f"{OUT}/cat/{wc.sample}.contig2classification.txt",
        db=f"{CAT_DB}/db" if CAT_DB else "",
        tax=f"{CAT_DB}/tax" if CAT_DB else "",
    shell:
        "CAT_pack contigs -c {input.contigs} -d {params.db} -t {params.tax} "
        "-o {params.prefix} -n {threads} --force && "
        "CAT_pack add_names -i {params.c2c} -o {output.named} -t {params.tax} --only_official"


rule classify_contigs:
    input:
        contigs=lambda wc: input_contigs(wc),
    output:
        report=f"{OUT}/contigs/{{sample}}.contigs.kreport",
        calls=f"{OUT}/contigs/{{sample}}.contigs.kraken",
    threads: THREADS
    params:
        cmd=lambda wc, output, threads: _classify_contigs_cmd(wc, output, threads),
    shell:
        "{params.cmd}"


# Horizontal coverage (breadth) per contig — fraction of each contig's bases covered by
# reads. At large reference scale this is a better "is it really present?" signal than depth.
rule contig_breadth:
    input:
        bam=f"{OUT}/binning/{{sample}}/aln.sorted.bam",
    output:
        cov=f"{OUT}/binning/{{sample}}/breadth.txt",
    shell:
        "samtools coverage {input.bam} -o {output.cov}"


# Score read-level kraken2 calls against the contig call (assembly as ground truth):
# concordant / discordant(misclassified) / unclassified, per contig and overall.
rule read_contig_accuracy:
    input:
        bam=f"{OUT}/binning/{{sample}}/aln.sorted.bam",
        read_kraken=f"{OUT}/kraken2/{{sample}}.{READ_LABEL}.kraken",
        contig_kraken=f"{OUT}/contigs/{{sample}}.contigs.kraken",
    output:
        tsv=f"{OUT}/reconcile/{{sample}}.read_accuracy.tsv",
        json=f"{OUT}/reconcile/{{sample}}.read_accuracy.json",
    params:
        sample=lambda wc: wc.sample,
        nodes=os.path.join(DB["kraken2"], "taxonomy", "nodes.dmp"),
    script:
        "../scripts/read_contig_accuracy.py"


rule reconcile:
    input:
        contig_calls=f"{OUT}/contigs/{{sample}}.contigs.kraken",
        depth=f"{OUT}/binning/{{sample}}/depth.txt",
        breadth=f"{OUT}/binning/{{sample}}/breadth.txt",
        read_report=f"{OUT}/kraken2/{{sample}}.{READ_LABEL}.kreport",
        cat_named=_cat_named,
    output:
        per_contig=f"{OUT}/reconcile/{{sample}}.contig_taxonomy.tsv",
        per_taxon=f"{OUT}/reconcile/{{sample}}.reconciliation.tsv",
        flags=f"{OUT}/reconcile/{{sample}}.flags.tsv",
        summary=f"{OUT}/reconcile/{{sample}}.reconcile.json",
    params:
        sample=lambda wc: wc.sample,
        read_label=READ_LABEL,
        cat_named=lambda wc: (f"{OUT}/cat/{wc.sample}.named.txt" if CAT_DB else ""),
    script:
        "../scripts/reconcile.py"
