import gzip
from metagx import formats, subsample, readfilter


def _write_fasta(p, n, length=100):
    with open(p, "w") as fh:
        for i in range(n):
            fh.write(f">read_{i}\n{'A'*length}\n")


def test_read_format(tmp_path):
    fa = tmp_path / "x.fasta"; fa.write_text(">a\nACGT\n")
    fq = tmp_path / "x.fastq"; fq.write_text("@a\nACGT\n+\nIIII\n")
    gz = tmp_path / "x.fq.gz"
    with gzip.open(gz, "wt") as fh:
        fh.write("@a\nACGT\n+\nIIII\n")
    assert formats.read_format(str(fa)) == "fasta"
    assert formats.read_format(str(fq)) == "fastq"
    assert formats.read_format(str(gz)) == "fastq"
    assert formats.is_gzipped(str(gz))


def test_estimate_read_length(tmp_path):
    fa = tmp_path / "r.fasta"; _write_fasta(str(fa), 50, length=120)
    est = formats.estimate_read_length(str(fa))
    assert est["median"] == 120 and est["n"] == 50


def test_subsample_seeded_reproducible(tmp_path):
    src = tmp_path / "in.fasta"; _write_fasta(str(src), 1000)
    o1, o2 = tmp_path / "o1.fasta", tmp_path / "o2.fasta"
    s1 = subsample.subsample(str(src), str(o1), 0.3, seed=42)
    s2 = subsample.subsample(str(src), str(o2), 0.3, seed=42)
    assert s1["kept"] == s2["kept"] and 0 < s1["kept"] < 1000
    assert o1.read_text() == o2.read_text()       # deterministic


def test_readfilter_include_exclude(tmp_path):
    # 3 reads: r1->taxid 10, r2->taxid 20, r3 unclassified(0)
    reads = tmp_path / "r.fasta"
    reads.write_text(">r1\nAAAA\n>r2\nCCCC\n>r3\nGGGG\n")
    kr = tmp_path / "r.kraken"
    kr.write_text("C\tr1\t10\t4\t\nC\tr2\t20\t4\t\nU\tr3\t0\t4\t\n")
    out = tmp_path / "f.fasta"
    # exclude taxid 10, keep unclassified -> keeps r2, r3
    st = readfilter.filter_reads(str(reads), str(kr), str(out), [10],
                                 mode="exclude", keep_unclassified=True)
    assert st["kept"] == 2 and st["removed"] == 1
    body = out.read_text()
    assert ">r1" not in body and ">r2" in body and ">r3" in body
    # include taxid 20 -> keeps only r2
    out2 = tmp_path / "f2.fasta"
    st2 = readfilter.filter_reads(str(reads), str(kr), str(out2), [20],
                                  mode="include", keep_unclassified=False)
    assert st2["kept"] == 1 and ">r2" in out2.read_text()
