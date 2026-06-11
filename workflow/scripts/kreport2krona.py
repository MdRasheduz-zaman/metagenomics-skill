"""Convert a kraken2 report into Krona text (for ktImportText).

A kraken2 report already encodes the full taxonomy: the name column is indented two spaces per
rank, so the indentation gives each taxon's lineage directly. We emit one row per taxon —
``<reads_at_taxon>\t<root>\t...\t<taxon>`` — which ``ktImportText`` aggregates up the tree.

This is what KrakenTools' kreport2krona.py does; we reimplement the ~15 lines so the workflow
needs no KrakenTools dependency, and — crucially — it works with **custom kraken2 DBs** whose
synthetic taxids are absent from any NCBI taxonomy (ktImportTaxonomy would drop them). Pure stdlib.
"""


def convert(kreport_path: str, out_path: str) -> int:
    """Write Krona text for one kraken2 report. Returns the number of taxa emitted."""
    lineage: list[str] = []
    n = 0
    with open(out_path, "w") as out:
        for line in open(kreport_path):
            c = line.rstrip("\n").split("\t")
            if len(c) < 6:
                continue
            taxon_reads, name = c[2], c[5]
            depth = (len(name) - len(name.lstrip(" "))) // 2   # 2 spaces per rank
            nm = name.strip()
            if nm.lower() == "unclassified":
                continue
            lineage = lineage[:depth]            # pop back to this taxon's parent
            lineage.append(nm)
            if taxon_reads.isdigit() and int(taxon_reads) > 0:
                out.write(taxon_reads + "\t" + "\t".join(lineage) + "\n")
                n += 1
    return n


if "snakemake" in globals():  # pragma: no cover - exercised inside the workflow
    convert(snakemake.input.kreport, snakemake.output.txt)  # noqa: F821
