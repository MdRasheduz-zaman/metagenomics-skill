# Optional random subsampling (seeded, format-aware, no external tools).
# Single-end only; paired-end is rejected in common.smk for now. Output keeps the input
# format (FASTA stays FASTA), uncompressed.

_FRACTION = (SUBSAMPLE or {}).get("fraction", 1.0)
_SEED = (SUBSAMPLE or {}).get("seed", 42)


rule subsample_fasta_se:
    wildcard_constraints:
        sample=_alt([s for s in FASTA_SAMPLES if layout(s) == "se"]),
    input:
        r1=lambda wc: source_reads(wc.sample)[0],
    output:
        r1=f"{OUT}/subsampled/{{sample}}_R1.fasta",
    run:
        stats = _subsample.subsample(input.r1, output.r1, _FRACTION, _SEED)
        print(f"[subsample] {wildcards.sample}: kept {stats['kept']}/{stats['total']} reads")


rule subsample_fastq_se:
    wildcard_constraints:
        sample=_alt([s for s in FASTQ_SAMPLES if layout(s) == "se"]),
    input:
        r1=lambda wc: source_reads(wc.sample)[0],
    output:
        r1=f"{OUT}/subsampled/{{sample}}_R1.fastq",
    run:
        stats = _subsample.subsample(input.r1, output.r1, _FRACTION, _SEED)
        print(f"[subsample] {wildcards.sample}: kept {stats['kept']}/{stats['total']} reads")
