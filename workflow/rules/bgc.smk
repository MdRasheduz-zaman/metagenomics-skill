# Biosynthetic gene cluster mining (Tier 2 functional add-on), gated by modules.bgc (requires
# assembly; WGS-only). antiSMASH detects secondary-metabolite clusters on assembled contigs.
# Isolated conda env (bgc.yaml); reference DBs via db.antismash (download-antismash-databases).
# antiSMASH writes a results tree under bgc/<sample>/; a sentinel marks completion.

rule antismash:
    wildcard_constraints:
        sample=_alt(WGS_SAMPLES),
    input:
        contigs=f"{OUT}/assembly/{{sample}}/final.contigs.fa",
    output:
        done=f"{OUT}/bgc/{{sample}}/antismash.done",
    threads: THREADS
    conda:
        "../envs/bgc.yaml"
    params:
        outdir=lambda wc: f"{OUT}/bgc/{wc.sample}",
        args=lambda wc, threads: " ".join(registry.render_args(
            "antismash", config.get("antismash", {}),
            managed={"cpus": threads, "output_dir": f"{OUT}/bgc/{wc.sample}",
                     "databases": DB.get("antismash", "")})),
    shell:
        "rm -rf {params.outdir} && antismash {params.args} {input.contigs} && "
        "touch {output.done}"
