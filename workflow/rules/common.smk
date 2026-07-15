# Shared setup: sample sheet parsing, the sweep axis, and registry-driven helpers.
# Imported by the top-level Snakefile before any rule modules.

import os
import re
from metagx import formats, registry
from metagx import subsample as _subsample

# --------------------------------------------------------------------------- #
# Paths / globals                                                             #
# --------------------------------------------------------------------------- #
PROJECT = config.get("project", "run")
OUTDIR = config.get("outdir", "results")
OUT = os.path.join(OUTDIR, PROJECT)
THREADS = int(config.get("threads", 8))
MODULES = config.get("modules", {})
DB = config.get("db", {})
# An optional db.build block builds the kraken2+Bracken DB as a pipeline step; classify/
# abundance then depend on the produced manifest so the build runs first (and only once).
DB_BUILD = DB.get("build")
# Mirror build_config: a db.build with no explicit db.kraken2 writes the DB to
# <outdir>/dbs/<strategy>. Default it here too so a hand-written/edited config that sets only
# db.build still runs, instead of KeyError'ing on DB["kraken2"] (RT-2).
if DB_BUILD and not DB.get("kraken2"):
    DB["kraken2"] = os.path.join(OUTDIR, "dbs", str(DB_BUILD.get("strategy", "standard")))
    DB.setdefault("bracken", DB["kraken2"])
DB_MANIFEST = os.path.join(DB["kraken2"], ".metagx_db.json") if DB.get("kraken2") else None


def db_ready_input():
    """The built-DB manifest, as a rule input — only when db.build is configured AND auto-build
    is on (db.build.auto, default True). With ``auto: false`` a plain run won't trigger a
    multi-hour build: the user builds it first (`metagx build-db` / target build_kraken2_db)."""
    auto = (DB_BUILD or {}).get("auto", True)
    return [DB_MANIFEST] if (DB_BUILD and DB_MANIFEST and auto) else []


# Per-tool module DBs (genomad/checkv/checkm2/gtdbtk/bakta/...) to auto-provision before use.
DB_PROVISION = list(DB.get("provision", []))


def provision_ready(tool):
    """A provision sentinel as a rule input — only when this tool is in db.provision (else [])
    so the domain/functional rule waits for its reference DB to be fetched (idempotently)."""
    return [os.path.join(OUT, "dbprovision", f"{tool}.done")] if tool in DB_PROVISION else []


# Target domains for the (optional) domain-taxonomy layer: viral / prokaryote / eukaryote.
DOMAINS = [d.lower() for d in config.get("domains", [])]


# --------------------------------------------------------------------------- #
# Samples: a TSV path (sample, r1, r2, platform, layout) or inline records.    #
#   platform: illumina|mgi (short) | ont (long) | pacbio_hifi|pacbio_clr (long) #
#   layout:   se | pe | interleaved  (default: pe if r2 present else se)        #
#   contigs:  optional pre-assembled FASTA (isolate genome / prior assembly /   #
#             references). When set, the assembler is skipped and these contigs #
#             feed the contig modules (functional/AMR, BGC, binning, domain).   #
#             A sample needs either reads (r1) or contigs.                       #
# --------------------------------------------------------------------------- #
SHORT_PLATFORMS = {"illumina", "mgi", "bgi"}
ONT_PLATFORMS = {"ont", "nanopore"}
PB_PLATFORMS = {"pacbio_hifi", "pacbio_clr", "pacbio"}
LONG_PLATFORMS = ONT_PLATFORMS | PB_PLATFORMS
KNOWN_PLATFORMS = SHORT_PLATFORMS | LONG_PLATFORMS


def _load_samples(spec):
    records = []
    if isinstance(spec, str):
        # BOM/whitespace-tolerant (Excel-saved sheets) — shared with the CLI/probe readers.
        records = formats.read_tsv_dicts(spec)
    else:
        records = spec
    samples = {}
    for rec in records:
        name = rec["sample"]
        r1 = rec.get("r1") or rec.get("R1")
        r2 = (rec.get("r2") or rec.get("R2") or "").strip() or None
        platform = (rec.get("platform") or "illumina").strip().lower()
        layout = (rec.get("layout") or "").strip().lower() or ("pe" if r2 else "se")
        library = (rec.get("library") or "wgs").strip().lower()
        # optional hybrid: long reads paired with a short-read sample (metaSPAdes only)
        long_reads = (rec.get("long_reads") or "").strip() or None
        long_platform = (rec.get("long_platform") or "ont").strip().lower()
        # optional: negative/blank control (decontam) and per-sample Bracken read length
        control = str(rec.get("control", "")).strip().lower() in {"1", "true", "yes", "y"}
        brl = (rec.get("bracken_read_length") or "").strip()
        bracken_read_length = int(brl) if brl else None
        # optional: group label for differential abundance (e.g. case/control)
        group = (rec.get("group") or "").strip() or None
        # optional: pre-assembled contigs/genome/MAGs (FASTA). When given, the assembler is
        # skipped and these contigs feed the assembly-dependent modules (functional/AMR, BGC,
        # binning, domain taxonomy, reconcile) directly — for users who already have an
        # isolate genome, a previous assembly, or downloaded references.
        contigs = (rec.get("contigs") or "").strip() or None
        samples[name] = {"r1": r1, "r2": r2, "platform": platform,
                         "layout": layout, "library": library,
                         "long_reads": long_reads, "long_platform": long_platform,
                         "control": control, "bracken_read_length": bracken_read_length,
                         "group": group, "contigs": contigs}
    if not samples:
        raise ValueError("No samples found.")
    return samples


SAMPLES = _load_samples(config["samples"])


def platform_of(sample):
    return SAMPLES[sample]["platform"]


def layout(sample):
    return SAMPLES[sample]["layout"]


def is_long(sample):
    return platform_of(sample) in LONG_PLATFORMS


def provided_contigs(sample):
    """User-supplied pre-assembled contigs/genome FASTA for this sample, or None."""
    return SAMPLES[sample].get("contigs")


def has_reads(sample):
    return bool(SAMPLES[sample].get("r1"))


def sample_format(sample):
    # A contigs-only sample (genome/MAGs provided, no reads) is FASTA by definition.
    if not has_reads(sample):
        return "fasta"
    return formats.read_format(SAMPLES[sample]["r1"])


def library_of(sample):
    return SAMPLES[sample]["library"]


def is_amplicon(sample):
    return library_of(sample) == "amplicon"


def is_ancient(sample):
    """Ancient/degraded DNA library — short-read shotgun with read-merging + damage auth."""
    return library_of(sample) == "ancient"


def is_control(sample):
    """Negative / blank control sample (used by the decontam module)."""
    return bool(SAMPLES[sample].get("control"))


def group_of(sample):
    """Group label for differential abundance (e.g. case/control), or None."""
    return SAMPLES[sample].get("group")


# --- validate the combinations up front so failures are clear, not cryptic ---
for _s in SAMPLES:
    if not has_reads(_s) and not provided_contigs(_s):
        raise ValueError(f"sample '{_s}': needs either reads (r1) or pre-assembled contigs.")
    if library_of(_s) not in {"wgs", "amplicon", "ancient"}:
        raise ValueError(f"sample '{_s}': library must be 'wgs', 'amplicon', or 'ancient'")
    if is_ancient(_s) and platform_of(_s) not in SHORT_PLATFORMS:
        raise ValueError(
            f"sample '{_s}': ancient library is short-read shotgun; platform must be short "
            f"(illumina/mgi/bgi), got '{platform_of(_s)}'"
        )
    if platform_of(_s) not in KNOWN_PLATFORMS:
        raise ValueError(
            f"sample '{_s}': unknown platform '{platform_of(_s)}'. "
            f"Known: {', '.join(sorted(KNOWN_PLATFORMS))}"
        )
    if layout(_s) not in {"se", "pe", "interleaved"}:
        raise ValueError(f"sample '{_s}': layout must be se|pe|interleaved")
    if is_long(_s) and (layout(_s) != "se" or SAMPLES[_s]["r2"]):
        raise ValueError(
            f"sample '{_s}': long-read platforms ({platform_of(_s)}) are single-end; "
            "drop r2 and use layout=se."
        )
    if layout(_s) == "interleaved" and platform_of(_s) not in SHORT_PLATFORMS:
        raise ValueError(f"sample '{_s}': interleaved layout applies to short reads only.")
    if layout(_s) == "pe" and not SAMPLES[_s]["r2"]:
        raise ValueError(f"sample '{_s}': layout=pe but no r2 given.")


def _grp(pred):
    return [s for s in SAMPLES if pred(s)]


FASTQ_SAMPLES = _grp(lambda s: sample_format(s) == "fastq")
FASTA_SAMPLES = _grp(lambda s: sample_format(s) == "fasta")
# Library strategy: assembly applies only to WGS (shotgun). Amplicon (marker-gene) reads
# are many copies of one locus, so assembly/binning/reconcile/etc. do not apply — they get
# the amplicon branch instead (primer trim -> OTU/abundance).
# WGS includes ancient (both are shotgun; ancient just adds read-merging + damage auth).
WGS_SAMPLES = _grp(lambda s: not is_amplicon(s))
AMPLICON_SAMPLES = _grp(is_amplicon)
ANCIENT_SAMPLES = _grp(is_ancient)
CONTROL_SAMPLES = _grp(is_control)

# Which QC tool cleans a sample: cutadapt (amplicon primers) | fastp_ancient (collapse
# overlapping ancient pairs) | fastp (short WGS) | ont | pacbio | None (FASTA WGS).
def qc_tool(sample):
    if is_amplicon(sample):
        return "cutadapt"
    if sample_format(sample) == "fasta":
        return None
    # ancient paired-end reads are merged (collapsed) into single reads first
    if is_ancient(sample) and layout(sample) == "pe":
        return "fastp_ancient"
    p = platform_of(sample)
    if p in SHORT_PLATFORMS:
        return "fastp"
    if p in ONT_PLATFORMS:
        return "ont"
    if p in PB_PLATFORMS:
        return "pacbio"
    return None


# Per-(tool, layout) sample sets used as wildcard constraints in qc.smk.
SHORT_PE_FASTQ = _grp(lambda s: qc_tool(s) == "fastp" and layout(s) == "pe")
SHORT_SE_FASTQ = _grp(lambda s: qc_tool(s) == "fastp" and layout(s) == "se")
SHORT_IL_FASTQ = _grp(lambda s: qc_tool(s) == "fastp" and layout(s) == "interleaved")
# Ancient reads merged into single collapsed reads (treated as SE downstream).
ANCIENT_MERGE = _grp(lambda s: qc_tool(s) == "fastp_ancient")
ONT_FASTQ = _grp(lambda s: qc_tool(s) == "ont")
PB_FASTQ = _grp(lambda s: qc_tool(s) == "pacbio")
# Amplicon QC = cutadapt primer removal (its own rule set, by layout).
AMPLICON_PE = _grp(lambda s: qc_tool(s) == "cutadapt" and layout(s) == "pe")
AMPLICON_SE = _grp(lambda s: qc_tool(s) == "cutadapt" and layout(s) != "pe")
# Amplicon profiling sets (post-QC): short -> VSEARCH OTUs or DADA2 ASVs, long -> Emu.
SHORT_AMPLICON = _grp(lambda s: is_amplicon(s) and platform_of(s) in SHORT_PLATFORMS)
LONG_AMPLICON = _grp(lambda s: is_amplicon(s) and platform_of(s) in LONG_PLATFORMS)
# Short-read amplicon denoising method: OTU clustering (VSEARCH, default) or ASVs (DADA2).
AMPLICON_METHOD = str(config.get("amplicon", {}).get("method", "otu")).lower()
SHORT_AMPLICON_OTU = SHORT_AMPLICON if AMPLICON_METHOD == "otu" else []
SHORT_AMPLICON_ASV = SHORT_AMPLICON if AMPLICON_METHOD == "asv" else []

# Sample -> group label map for differential abundance (only labelled samples).
GROUPS = {s: group_of(s) for s in SAMPLES if group_of(s)}

# Optional random subsampling (single-end only for now).
SUBSAMPLE = config.get("subsample")
if SUBSAMPLE:
    _bad = _grp(lambda s: layout(s) != "se")
    if _bad:
        raise ValueError(
            "subsample currently supports single-end only; offending samples: "
            f"{', '.join(_bad)} (paired/interleaved)."
        )


def _alt(names):
    """Regex alternation over names; a never-match pattern when the list is empty."""
    return "|".join(re.escape(n) for n in names) if names else r"(?!)"


def cat_cmd(path):
    """Stream a (possibly gzipped) file to stdout, portably."""
    return f"gunzip -c {path}" if str(path).endswith(".gz") else f"cat {path}"


# --------------------------------------------------------------------------- #
# Sweep axis: always present. From config.sweep, else a single point from the  #
# kraken2 section (or kraken2's own default).                                  #
# --------------------------------------------------------------------------- #
def _sweep():
    sw = config.get("sweep")
    if sw:
        return sw["param"], [str(v) for v in sw["values"]]
    param = "confidence"
    pinned = config.get("kraken2", {}).get(param)
    if pinned is None:
        pinned = registry.load_registry("kraken2")["params"][param].get("default", 0.0)
    return param, [str(pinned)]


SWEEP_PARAM, SWEEP_VALUES = _sweep()


def sweep_label(value):
    return f"{SWEEP_PARAM}_{value}"


SWEEP_LABELS = [sweep_label(v) for v in SWEEP_VALUES]
LABEL_TO_VALUE = dict(zip(SWEEP_LABELS, SWEEP_VALUES))
# Read-level confidence used as the reference when reconciling against contigs
# (the first / most-permissive sweep point).
READ_LABEL = SWEEP_LABELS[0]


# --------------------------------------------------------------------------- #
# Read routing chain:  raw -> (subsample) -> (QC by platform) -> classify       #
# --------------------------------------------------------------------------- #
def is_paired(sample):
    """True when classification should pass two read files (pe, or split interleaved).

    Ancient samples are collapsed to a single merged stream, so they are single-end
    downstream regardless of their input layout.
    """
    if qc_tool(sample) == "fastp_ancient":
        return False
    return layout(sample) in {"pe", "interleaved"}


def two_file_qc(sample):
    """QC emits two files for pe and (split) interleaved short reads; one otherwise.

    Ancient merge collapses pairs into one file, so it emits a single file.
    """
    if qc_tool(sample) == "fastp_ancient":
        return False
    return layout(sample) in {"pe", "interleaved"}


# --------------------------------------------------------------------------- #
# Per-sample Bracken read length: sample field > per-platform map > global.    #
# (Bracken's kmer_distrib is length-specific, so mixed short+long runs need    #
#  the right length per sample — the DB must have that databaseXmers built.)   #
# --------------------------------------------------------------------------- #
BRACKEN_CFG = config.get("bracken", {})
BRACKEN_LEN_BY_PLATFORM = config.get("bracken_read_length_by_platform", {})


def bracken_read_length(sample):
    s = SAMPLES[sample]
    if s.get("bracken_read_length"):
        return int(s["bracken_read_length"])
    p = platform_of(sample)
    if p in BRACKEN_LEN_BY_PLATFORM:
        return int(BRACKEN_LEN_BY_PLATFORM[p])
    return int(BRACKEN_CFG.get("read_length", 100))


def raw_reads(sample):
    s = SAMPLES[sample]
    return [s["r1"]] + ([s["r2"]] if s["r2"] else [])


# Optional host removal as a first-class pre-classification step (map to host, keep unmapped).
HOST_GENOME = config.get("host_removal", {}).get("genome")


def source_reads(sample):
    """Reads entering the pipeline: host-depleted (if host_removal.genome set, SE) else raw."""
    if HOST_GENOME and layout(sample) == "se":
        return [f"{OUT}/hostclean/{sample}.fastq"]
    return raw_reads(sample)


def qc_active(sample):
    return MODULES.get("qc", True) and qc_tool(sample) is not None


def subsample_out(sample):
    ext = formats.canonical_ext(SAMPLES[sample]["r1"])
    return [f"{OUT}/subsampled/{sample}_R1{ext}"]  # single-end only (validated above)


def staged_reads(sample):
    """Reads after host removal + optional subsampling — input to QC / classification."""
    return subsample_out(sample) if SUBSAMPLE else source_reads(sample)


def qc_out(sample):
    out = [f"{OUT}/qc/{sample}_R1.fastq.gz"]
    if two_file_qc(sample):
        out.append(f"{OUT}/qc/{sample}_R2.fastq.gz")
    return out


def reads_for_classify(sample):
    if qc_active(sample):
        return qc_out(sample)
    return staged_reads(sample)


# --------------------------------------------------------------------------- #
# Assembly tool + mapping preset chosen by platform                            #
# --------------------------------------------------------------------------- #
# Short-read assembler is configurable: megahit (default, fast) or metaspades
# (accurate, paired-end only, supports hybrid with long reads). Long-read samples
# always use Flye/metaFlye.
ASSEMBLY = config.get("assembly", {})
SHORT_ASSEMBLER = str(ASSEMBLY.get("assembler", "megahit")).lower()
if SHORT_ASSEMBLER not in {"megahit", "metaspades"}:
    raise ValueError(f"assembly.assembler must be megahit|metaspades, got '{SHORT_ASSEMBLER}'")


def long_reads_of(sample):
    return SAMPLES[sample].get("long_reads")


def is_hybrid(sample):
    """Short-read sample that also carries long reads for hybrid metaSPAdes assembly."""
    return bool(long_reads_of(sample)) and not is_long(sample)


def assembler(sample):
    if is_long(sample):
        return "flye"
    return SHORT_ASSEMBLER  # megahit | metaspades


def flye_read_flag(sample):
    return {"ont": "--nano-hq", "nanopore": "--nano-hq",
            "pacbio_hifi": "--pacbio-hifi", "pacbio": "--pacbio-hifi",
            "pacbio_clr": "--pacbio-raw"}.get(platform_of(sample), "--nano-hq")


def minimap2_preset(sample):
    p = platform_of(sample)
    if p in ONT_PLATFORMS:
        return "map-ont"
    if p == "pacbio_hifi":
        return "map-hifi"
    if p in PB_PLATFORMS:
        return "map-pb"
    return "sr"  # short reads


# Samples that supply their own contigs skip the assembler entirely — the staging rule
# (assembly.smk) puts their FASTA where the assembler's output would be.
PROVIDED_CONTIGS = _grp(lambda s: provided_contigs(s))

# Assembly is WGS-only (amplicon excluded — assembly does not apply to marker-gene data).
# Samples with provided contigs are excluded from every assembler so they don't also try to
# (re)assemble — the staged contigs satisfy the same output.
SHORT_ASM = _grp(lambda s: assembler(s) == "megahit" and not is_amplicon(s)
                 and s not in PROVIDED_CONTIGS)
SPADES_ASM = _grp(lambda s: assembler(s) == "metaspades" and not is_amplicon(s)
                  and s not in PROVIDED_CONTIGS)
LONG_ASM = _grp(lambda s: assembler(s) == "flye" and not is_amplicon(s)
                and s not in PROVIDED_CONTIGS)

# metaSPAdes needs paired-end short reads (single-end is unsupported, even in hybrid mode).
_spades_se = [s for s in SPADES_ASM if layout(s) == "se"]
if _spades_se:
    raise ValueError(
        "assembly.assembler=metaspades requires paired-end short reads; single-end "
        f"sample(s) {', '.join(_spades_se)} cannot use metaSPAdes (use megahit, or provide r2)."
    )


# Constrain wildcards so sample names and sweep labels are matched unambiguously.
wildcard_constraints:
    sample=_alt(list(SAMPLES)),
    label=_alt(SWEEP_LABELS),
    level="[A-Z][0-9]?",
