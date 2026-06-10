import pytest
from metagx import config_builder as cb
from metagx import registry

DB = {"kraken2": "k"}
S = [{"sample": "a", "r1": "a.fastq.gz"}]


def test_minimal_valid():
    cfg = cb.build_config(samples=S, db=DB, modules={"assembly": False})
    assert cfg["db"]["bracken"] == "k"   # defaults to kraken2
    assert cfg["modules"]["classify"] is True


def test_module_dependencies():
    with pytest.raises(registry.ValidationError):
        cb.build_config(samples=S, db=DB, modules={"binning": True, "assembly": False})
    with pytest.raises(registry.ValidationError):
        cb.build_config(samples=S, db=DB, modules={"stats": True, "abundance": False})
    with pytest.raises(registry.ValidationError):
        cb.build_config(samples=S, db=DB,
                        modules={"reconcile": True, "assembly": False, "classify": True})


def test_preset_merge_and_override():
    cfg = cb.build_config(samples=S, db=DB, preset="pathogen-detection",
                          bracken={"threshold": 50})
    assert cfg["bracken"]["threshold"] == 50          # user override wins
    assert cfg["bracken"]["read_length"] == 150       # from preset


def test_sweep_conflict():
    with pytest.raises(registry.ValidationError):
        cb.build_config(samples=S, db=DB, sweep={"param": "confidence", "values": [0.1]},
                        kraken2={"confidence": 0.2})


def test_all_amplicon_blocks_assembly():
    amp = [{"sample": "a", "r1": "a.fq", "library": "amplicon"}]
    with pytest.raises(registry.ValidationError):
        cb.build_config(samples=amp, db=DB, modules={"assembly": True})
    # amplicon-only without assembly is fine
    cfg = cb.build_config(samples=amp, db=DB, modules={"assembly": False})
    assert cfg is not None


def test_library_and_platform_validation():
    with pytest.raises(registry.ValidationError):
        cb.build_config(samples=[{"sample": "a", "r1": "x", "library": "bogus"}], db=DB)
    with pytest.raises(registry.ValidationError):
        cb.build_config(samples=[{"sample": "a", "r1": "x", "platform": "ont", "r2": "y"}], db=DB)


def test_subsample_validation():
    with pytest.raises(registry.ValidationError):
        cb.build_config(samples=S, db=DB, subsample={"fraction": 2.0})


def test_functional_module_and_db_extras():
    cfg = cb.build_config(
        samples=S,
        db={"kraken2": "k", "humann_nucleotide": "n", "amrfinderplus": "amr"},
        modules={"functional": True},
        abricate={"db": "card"},
    )
    assert cfg["modules"]["functional"] is True
    assert cfg["abricate"]["db"] == "card"
    assert cfg["db"]["humann_nucleotide"] == "n"   # functional db extras carried through
    assert cfg["db"]["amrfinderplus"] == "amr"


def test_functional_blocked_for_all_amplicon():
    amp = [{"sample": "a", "r1": "a.fq", "library": "amplicon"}]
    with pytest.raises(registry.ValidationError):
        cb.build_config(samples=amp, db=DB, modules={"functional": True})


def test_bin_refinement_dependencies():
    # refinement needs binning (which needs assembly)
    with pytest.raises(registry.ValidationError):
        cb.build_config(samples=S, db=DB,
                        modules={"bin_refinement": True, "binning": False})
    cfg = cb.build_config(
        samples=S, db=DB,
        modules={"assembly": True, "binning": True, "bin_refinement": True},
        das_tool={"score_threshold": 0.6}, drep={"s_ani": 0.99},
    )
    assert cfg["modules"]["bin_refinement"] is True
    assert cfg["das_tool"]["score_threshold"] == 0.6
    assert cfg["drep"]["s_ani"] == 0.99


SPE = [{"sample": "a", "r1": "a_1.fq.gz", "r2": "a_2.fq.gz"}]


def test_assembler_choice_metaspades():
    cfg = cb.build_config(samples=SPE, db=DB, modules={"assembly": True},
                          assembly={"assembler": "metaspades"}, metaspades={"memory_gb": 64})
    assert cfg["assembly"]["assembler"] == "metaspades"
    assert cfg["metaspades"]["memory_gb"] == 64
    with pytest.raises(registry.ValidationError):
        cb.build_config(samples=SPE, db=DB, assembly={"assembler": "bogus"})


def test_hybrid_long_reads_requires_metaspades():
    hyb = [{"sample": "a", "r1": "a_1.fq.gz", "r2": "a_2.fq.gz",
            "long_reads": "a_ont.fq.gz", "long_platform": "ont"}]
    with pytest.raises(registry.ValidationError):
        cb.build_config(samples=hyb, db=DB, modules={"assembly": True})  # default megahit
    cfg = cb.build_config(samples=hyb, db=DB, modules={"assembly": True},
                          assembly={"assembler": "metaspades"})
    assert cfg["assembly"]["assembler"] == "metaspades"
    # a long-read sample cannot also carry long_reads (it is already long)
    with pytest.raises(registry.ValidationError):
        cb.build_config(samples=[{"sample": "a", "r1": "a.fq", "platform": "ont",
                                  "long_reads": "x.fq"}], db=DB,
                        assembly={"assembler": "metaspades"})


def test_classify_consensus_dependencies():
    with pytest.raises(registry.ValidationError):
        cb.build_config(samples=S, db=DB,
                        modules={"classify_consensus": True, "classify": False})
    cfg = cb.build_config(samples=S, db={"kraken2": "k", "kaiju": "kdb"},
                          modules={"classify_consensus": True},
                          consensus={"classifier": "kaiju"}, kaiju={"min_match_length": 12})
    assert cfg["consensus"]["classifier"] == "kaiju"
    assert cfg["kaiju"]["min_match_length"] == 12
    assert cfg["db"]["kaiju"] == "kdb"
    with pytest.raises(registry.ValidationError):
        cb.build_config(samples=S, db=DB, consensus={"classifier": "bogus"})


def test_aggregate_requires_classify():
    with pytest.raises(registry.ValidationError):
        cb.build_config(samples=S, db=DB, modules={"aggregate": True, "classify": False})
    cfg = cb.build_config(samples=S, db=DB, modules={"aggregate": True})
    assert cfg["modules"]["aggregate"] is True


ANC = [{"sample": "a", "r1": "a_1.fq.gz", "r2": "a_2.fq.gz", "library": "ancient"}]


def test_ancient_library_and_damage():
    # damage needs assembly + an ancient sample
    with pytest.raises(registry.ValidationError):
        cb.build_config(samples=ANC, db=DB, modules={"damage": True, "assembly": False})
    with pytest.raises(registry.ValidationError):
        cb.build_config(samples=S, db=DB, modules={"damage": True, "assembly": True})  # no ancient
    cfg = cb.build_config(samples=ANC, db=DB, modules={"assembly": True, "damage": True},
                          mapdamage={"length": 50})
    assert cfg["modules"]["damage"] is True
    assert cfg["mapdamage"]["length"] == 50
    # ancient must be short-read
    with pytest.raises(registry.ValidationError):
        cb.build_config(samples=[{"sample": "a", "r1": "a.fq", "platform": "ont",
                                  "library": "ancient"}], db=DB)


def test_strain_requires_assembly():
    with pytest.raises(registry.ValidationError):
        cb.build_config(samples=S, db=DB, modules={"strain": True, "assembly": False})
    cfg = cb.build_config(samples=S, db=DB, modules={"assembly": True, "strain": True},
                          instrain={"min_cov": 10})
    assert cfg["instrain"]["min_cov"] == 10


def test_decontam_needs_abundance_and_control():
    with pytest.raises(registry.ValidationError):
        cb.build_config(samples=S, db=DB, modules={"decontam": True, "abundance": False})
    with pytest.raises(registry.ValidationError):
        cb.build_config(samples=S, db=DB, modules={"decontam": True})  # no control sample
    ctrls = [{"sample": "blank", "r1": "b.fq.gz", "control": True},
             {"sample": "s1", "r1": "s1.fq.gz"}]
    cfg = cb.build_config(samples=ctrls, db=DB, modules={"decontam": True})
    assert cfg["modules"]["decontam"] is True


def test_per_sample_bracken_length():
    sheet = [{"sample": "s1", "r1": "s1.fq.gz", "bracken_read_length": 150}]
    cfg = cb.build_config(samples=sheet, db=DB,
                          bracken_read_length_by_platform={"ont": 1000})
    assert cfg["samples"][0]["bracken_read_length"] == 150
    assert cfg["bracken_read_length_by_platform"]["ont"] == 1000
    with pytest.raises(registry.ValidationError):
        cb.build_config(samples=S, db=DB, bracken_read_length_by_platform={"bogus": 100})


def test_ancient_dna_preset():
    cfg = cb.build_config(samples=ANC, db=DB, preset="ancient-dna")
    assert cfg["modules"]["damage"] is True
    assert cfg["fastp"]["merge"] is True       # preset enables read-merging
    assert cfg["bracken"]["read_length"] == 50


GRP = [{"sample": "a1", "r1": "a1.fq", "group": "case"},
       {"sample": "a2", "r1": "a2.fq", "group": "case"},
       {"sample": "b1", "r1": "b1.fq", "group": "control"},
       {"sample": "b2", "r1": "b2.fq", "group": "control"}]


def test_differential_needs_abundance_and_two_groups():
    with pytest.raises(registry.ValidationError):
        cb.build_config(samples=GRP, db=DB, modules={"differential": True, "abundance": False})
    # one group only -> rejected
    one = [{"sample": "a", "r1": "a.fq", "group": "case"},
           {"sample": "b", "r1": "b.fq", "group": "case"}]
    with pytest.raises(registry.ValidationError):
        cb.build_config(samples=one, db=DB, modules={"differential": True})
    # under-replicated group -> rejected (permutation needs >=2 per group)
    thin = [{"sample": "a", "r1": "a.fq", "group": "case"},
            {"sample": "b", "r1": "b.fq", "group": "control"}]
    with pytest.raises(registry.ValidationError):
        cb.build_config(samples=thin, db=DB, modules={"differential": True})
    cfg = cb.build_config(samples=GRP, db=DB, modules={"differential": True},
                          differential={"n_permutations": 499, "fdr": 0.1})
    assert cfg["modules"]["differential"] is True
    assert cfg["differential"]["n_permutations"] == 499
    assert cfg["differential"]["fdr"] == 0.1


def test_differential_tuning_validation():
    with pytest.raises(registry.ValidationError):
        cb.build_config(samples=GRP, db=DB, modules={"differential": True},
                        differential={"n_permutations": 10})       # too few
    with pytest.raises(registry.ValidationError):
        cb.build_config(samples=GRP, db=DB, modules={"differential": True},
                        differential={"fdr": 1.5})                 # out of range


def test_bgc_requires_assembly_and_is_wgs_only():
    with pytest.raises(registry.ValidationError):
        cb.build_config(samples=SPE, db=DB, modules={"bgc": True, "assembly": False})
    amp = [{"sample": "a", "r1": "a.fq", "library": "amplicon"}]
    with pytest.raises(registry.ValidationError):
        cb.build_config(samples=amp, db=DB, modules={"bgc": True, "assembly": True})
    cfg = cb.build_config(samples=SPE, db={"kraken2": "k", "antismash": "as"},
                          modules={"assembly": True, "bgc": True},
                          antismash={"taxon": "bacteria", "cb_general": True})
    assert cfg["modules"]["bgc"] is True
    assert cfg["antismash"]["taxon"] == "bacteria"
    assert cfg["db"]["antismash"] == "as"           # antismash db carried through
    with pytest.raises(registry.ValidationError):
        cb.build_config(samples=SPE, db=DB, modules={"assembly": True, "bgc": True},
                        antismash={"taxon": "martian"})            # bad enum


def test_amplicon_method_asv_and_dada2():
    amp = [{"sample": "x", "r1": "x_1.fq", "r2": "x_2.fq", "library": "amplicon"}]
    cfg = cb.build_config(samples=amp, db=DB, modules={"assembly": False},
                          amplicon={"fwd_primer": "AAA", "rev_primer": "TTT", "method": "asv"},
                          dada2={"trunc_len_f": 240, "max_ee_f": 2.0})
    assert cfg["amplicon"]["method"] == "asv"
    assert cfg["dada2"]["trunc_len_f"] == 240
    # default method is otu
    cfg2 = cb.build_config(samples=amp, db=DB, modules={"assembly": False},
                           amplicon={"fwd_primer": "AAA", "rev_primer": "TTT"})
    assert cfg2["amplicon"]["method"] == "otu"
    with pytest.raises(registry.ValidationError):
        cb.build_config(samples=amp, db=DB, modules={"assembly": False},
                        amplicon={"method": "bogus"})


def test_amr_preset_merges_new_tool_sections():
    # Generic preset merge must carry tool sections for tools beyond the original eight.
    cfg = cb.build_config(samples=S, db=DB, preset="amr-surveillance",
                          abricate={"minid": 95})
    assert cfg["modules"]["functional"] is True
    assert cfg["abricate"]["db"] == "card"        # from preset
    assert cfg["abricate"]["minid"] == 95         # user override wins
    assert cfg["amrfinderplus"]["plus"] is True   # preset-only tool section preserved
