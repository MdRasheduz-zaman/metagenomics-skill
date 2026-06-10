# Bracken abundance re-estimation, per (sample, sweep value), at the configured level.

BRACKEN_LEVEL = config.get("bracken", {}).get("level", "S")


def _bracken_cmd(wc, input, output):
    base = {k: v for k, v in config.get("bracken", {}).items()}
    base.pop("read_length_by_platform", None)  # not a Bracken flag; routing hint only
    # per-sample read length (short mate length / long median); DB must have its databaseXmers
    base["read_length"] = bracken_read_length(wc.sample)
    managed = {
        "db": DB.get("bracken") or DB["kraken2"],
        "input": input.report,
        "output": output.bracken,
        "report_out": output.report,
    }
    args = registry.render_args("bracken", base, managed=managed)
    return "bracken " + " ".join(args)


rule bracken:
    input:
        report=f"{OUT}/kraken2/{{sample}}.{{label}}.kreport",
    output:
        bracken=f"{OUT}/bracken/{{sample}}.{{label}}.{BRACKEN_LEVEL}.bracken",
        report=f"{OUT}/bracken/{{sample}}.{{label}}.{BRACKEN_LEVEL}.breport",
    params:
        cmd=lambda wc, input, output: _bracken_cmd(wc, input, output),
    shell:
        "{params.cmd}"


rule bracken_combined:
    input:
        tables=expand(
            f"{OUT}/bracken/{{sample}}.{{label}}.{BRACKEN_LEVEL}.bracken",
            sample=list(SAMPLES), label=SWEEP_LABELS,
        ),
    output:
        tsv=f"{OUT}/summary/bracken_combined.tsv",
    script:
        "../scripts/combine_bracken.py"
