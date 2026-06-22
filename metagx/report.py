"""Provenance capture + Methods/report generation ("Copy as Methods").

Given a config and a finished (or partial) results directory, this builds:
  * a machine-readable run manifest (tool versions, exact commands, db identity,
    QC + classification metrics, timestamps)
  * a human Methods paragraph with citations
  * a full Markdown report (methods + parameter table + figures + abundance table),
    optionally rendered to LaTeX / PDF.

Everything is derived from the parameter registries, so the reported commands are exactly
what the workflow runs.
"""

from __future__ import annotations

import csv
import datetime as _dt
import glob
import hashlib
import json
import os
import re
import shlex
import shutil
import subprocess
from typing import Any, Dict, List, Optional

import yaml

from . import __version__, registry

# Citations keyed by tool. Kept here so the Methods section is publication-ready.
CITATIONS: Dict[str, str] = {
    "fastp": "Chen S, et al. fastp: an ultra-fast all-in-one FASTQ preprocessor. "
             "Bioinformatics. 2018;34(17):i884-i890.",
    "kraken2": "Wood DE, Lu J, Langmead B. Improved metagenomic analysis with Kraken 2. "
               "Genome Biology. 2019;20:257.",
    "bracken": "Lu J, et al. Bracken: estimating species abundance in metagenomics data. "
               "PeerJ Computer Science. 2017;3:e104.",
    "megahit": "Li D, et al. MEGAHIT: an ultra-fast single-node solution for large and "
               "complex metagenomics assembly. Bioinformatics. 2015;31(10):1674-1676.",
    "metaspades": "Nurk S, et al. metaSPAdes: a new versatile metagenomic assembler. "
                  "Genome Research. 2017;27(5):824-834.",
    "metabat2": "Kang DD, et al. MetaBAT 2: an adaptive binning algorithm for robust and "
                "efficient genome reconstruction from metagenome assemblies. PeerJ. 2019;7:e7359.",
    "minimap2": "Li H. Minimap2: pairwise alignment for nucleotide sequences. "
                "Bioinformatics. 2018;34(18):3094-3100.",
    "samtools": "Danecek P, et al. Twelve years of SAMtools and BCFtools. GigaScience. 2021;10(2).",
    "snakemake": "Mölder F, et al. Sustainable data analysis with Snakemake. "
                 "F1000Research. 2021;10:33.",
    "porechop_abi": "Bonenfant Q, et al. Porechop_ABI: discovering unknown adapters in ONT "
                    "sequencing reads. Bioinformatics Advances. 2023;3(1):vbac085.",
    "chopper": "De Coster W, Rademakers R. NanoPack2: population-scale evaluation of "
               "long-read sequencing data. Bioinformatics. 2023;39(5):btad311.",
    "flye": "Kolmogorov M, et al. metaFlye: scalable long-read metagenome assembly using "
            "repeat graphs. Nature Methods. 2020;17:1103-1110.",
    "genomad": "Camargo AP, et al. Identification of mobile genetic elements with geNomad. "
               "Nature Biotechnology. 2024;42:1303-1312.",
    "checkv": "Nayfach S, et al. CheckV assesses the quality and completeness of "
              "metagenome-assembled viral genomes. Nature Biotechnology. 2021;39:578-585.",
    "checkm2": "Chklovski A, et al. CheckM2: a rapid, scalable and accurate tool for "
               "assessing microbial genome quality using machine learning. Nature Methods. "
               "2023;20:1203-1212.",
    "gtdbtk": "Chaumeil PA, et al. GTDB-Tk v2: memory friendly classification with the "
              "Genome Taxonomy Database. Bioinformatics. 2022;38(23):5315-5316.",
    "eukrep": "West PT, et al. Genome-reconstruction for eukaryotes from complex natural "
              "microbial communities. Genome Research. 2018;28:569-580.",
    "eukcc": "Saary P, et al. Estimating the quality of eukaryotic genomes recovered from "
             "metagenomic analysis with EukCC. Genome Biology. 2020;21:244.",
    "humann": "Beghini F, et al. Integrating taxonomic, functional, and strain-level profiling "
              "of diverse microbial communities with bioBakery 3. eLife. 2021;10:e65088.",
    "metaphlan": "Blanco-Míguez A, et al. Extending and improving metagenomic taxonomic "
                 "profiling with uncharacterized species using MetaPhlAn 4. Nature "
                 "Biotechnology. 2023;41:1633-1644.",
    "amrfinderplus": "Feldgarden M, et al. AMRFinderPlus and the Reference Gene Catalog "
                     "facilitate examination of the genomic links among antimicrobial "
                     "resistance, stress response, and virulence. Scientific Reports. 2021;11:12728.",
    "abricate": "Seemann T. ABRicate: mass screening of contigs for antimicrobial and "
                "virulence genes. https://github.com/tseemann/abricate.",
    "bakta": "Schwengers O, et al. Bakta: rapid and standardized annotation of bacterial "
             "genomes via alignment-free sequence identification. Microbial Genomics. 2021;7(11).",
    "eggnog": "Cantalapiedra CP, et al. eggNOG-mapper v2: functional annotation, orthology "
              "assignments, and domain prediction at the metagenomic scale. Molecular Biology "
              "and Evolution. 2021;38(12):5825-5829.",
    "maxbin2": "Wu YW, Simmons BA, Singer SW. MaxBin 2.0: an automated binning algorithm to "
               "recover genomes from multiple metagenomic datasets. Bioinformatics. "
               "2016;32(4):605-607.",
    "concoct": "Alneberg J, et al. Binning metagenomic contigs by coverage and composition. "
               "Nature Methods. 2014;11:1144-1146.",
    "das_tool": "Sieber CMK, et al. Recovery of genomes from metagenomes via a dereplication, "
                "aggregation and scoring strategy. Nature Microbiology. 2018;3:836-843.",
    "drep": "Olm MR, et al. dRep: a tool for fast and accurate genomic comparisons that "
            "enables improved genome recovery from metagenomes. ISME J. 2017;11:2864-2868.",
    "kaiju": "Menzel P, Ng KL, Krogh A. Fast and sensitive taxonomic classification for "
             "metagenomics with Kaiju. Nature Communications. 2016;7:11257.",
    "multiqc": "Ewels P, et al. MultiQC: summarize analysis results for multiple tools and "
               "samples in a single report. Bioinformatics. 2016;32(19):3047-3048.",
    "krona": "Ondov BD, Bergman NH, Phillippy AM. Interactive metagenomic visualization in a "
             "Web browser. BMC Bioinformatics. 2011;12:385.",
    "mapdamage": "Jónsson H, et al. mapDamage2.0: fast approximate Bayesian estimates of "
                 "ancient DNA damage parameters. Bioinformatics. 2013;29(13):1682-1684.",
    "instrain": "Olm MR, et al. inStrain profiles population microdiversity from metagenomic "
                "data and sensitively detects shared microbial strains. Nature Biotechnology. "
                "2021;39:727-736.",
    "antismash": "Blin K, et al. antiSMASH 7.0: new and improved predictions for detection, "
                 "regulation, chemical structures and visualisation. Nucleic Acids Research. "
                 "2023;51(W1):W46-W50.",
    "dada2": "Callahan BJ, et al. DADA2: high-resolution sample inference from Illumina "
             "amplicon data. Nature Methods. 2016;13:581-583.",
    "vsearch": "Rognes T, et al. VSEARCH: a versatile open source tool for metagenomics. "
               "PeerJ. 2016;4:e2584.",
    "emu": "Curry KD, et al. Emu: species-level microbial community profiling of full-length "
           "16S rRNA Oxford Nanopore sequencing data. Nature Methods. 2022;19:845-853.",
    "cutadapt": "Martin M. Cutadapt removes adapter sequences from high-throughput "
                "sequencing reads. EMBnet.journal. 2011;17(1):10-12.",
}

MODULE_TOOLS = {
    "qc": ["fastp"],
    "classify": ["kraken2"],
    "classify_consensus": ["metaphlan", "kaiju"],
    "abundance": ["bracken"],
    "assembly": ["megahit"],
    "binning": ["minimap2", "samtools", "metabat2"],
    "bin_refinement": ["maxbin2", "concoct", "das_tool", "drep"],
    "functional": ["humann", "metaphlan", "amrfinderplus", "abricate", "bakta", "eggnog"],
    "aggregate": ["multiqc", "krona"],
    "damage": ["mapdamage"],
    "strain": ["instrain"],
    "bgc": ["antismash"],
}


# --------------------------------------------------------------------------- #
# Helpers                                                                      #
# --------------------------------------------------------------------------- #
def _outdir(cfg: Dict[str, Any]) -> str:
    return os.path.join(cfg.get("outdir", "results"), cfg.get("project", "run"))


SHORT_P = {"illumina", "mgi", "bgi"}
ONT_P = {"ont", "nanopore"}
PB_P = {"pacbio_hifi", "pacbio_clr", "pacbio"}


def _records(cfg: Dict[str, Any]) -> List[Dict[str, Any]]:
    """The run's sample records, whether given inline or as a TSV path."""
    samples = cfg.get("samples")
    if isinstance(samples, str) and os.path.isfile(samples):
        with open(samples) as fh:
            return list(csv.DictReader(fh, delimiter="\t"))
    if isinstance(samples, list):
        return samples
    return []


def _platforms(cfg: Dict[str, Any]) -> set:
    """Set of sequencing platforms across the run's samples (for tool selection)."""
    plats = {str(r.get("platform", "illumina")).lower() for r in _records(cfg)}
    return plats or {"illumina"}


def active_tools(cfg: Dict[str, Any]) -> List[str]:
    """Tools actually used, given enabled modules AND sample platforms."""
    mods = cfg.get("modules", {})
    plats = _platforms(cfg)
    tools: List[str] = []
    if mods.get("qc"):
        if plats & SHORT_P:
            tools.append("fastp")
        if plats & ONT_P:
            tools += ["porechop_abi", "chopper"]
        if plats & PB_P:
            tools.append("chopper")
    if mods.get("classify"):
        tools.append("kraken2")
    if mods.get("classify_consensus"):
        tools.append(str(cfg.get("consensus", {}).get("classifier", "metaphlan")).lower())
    if mods.get("abundance"):
        tools.append("bracken")
    if mods.get("assembly"):
        asm = str(cfg.get("assembly", {}).get("assembler", "megahit")).lower()
        if plats & SHORT_P:
            tools.append("metaspades" if asm == "metaspades" else "megahit")
        if plats & (ONT_P | PB_P):
            tools.append("flye")
    if mods.get("binning"):
        tools += ["minimap2", "samtools", "metabat2"]
    if mods.get("bin_refinement"):
        tools += ["maxbin2", "concoct", "das_tool", "drep"]
    if mods.get("domain_taxonomy"):
        doms = [d.lower() for d in cfg.get("domains", [])]
        if "viral" in doms:
            tools += ["genomad", "checkv"]
        if "prokaryote" in doms:
            tools += ["gtdbtk", "checkm2"]
        if "eukaryote" in doms:
            tools += ["eukrep", "eukcc"]
    if mods.get("functional"):
        # read-based pathway profiling always runs; AMR needs contigs, annotation needs bins
        tools += ["humann", "metaphlan"]
        if mods.get("assembly"):
            tools += ["amrfinderplus", "abricate"]
        if mods.get("binning"):
            tools += ["bakta", "eggnog"]
    if mods.get("aggregate"):
        tools += ["multiqc", "krona"]
    if mods.get("strain"):
        tools.append("instrain")
    if mods.get("damage"):
        tools.append("mapdamage")
    if mods.get("bgc") and mods.get("assembly"):
        tools.append("antismash")
    # amplicon branch is data-driven (per-sample library=amplicon), not a module toggle
    amp = [r for r in _records(cfg) if str(r.get("library", "wgs")).lower() == "amplicon"]
    if amp:
        tools.append("cutadapt")
        method = str(cfg.get("amplicon", {}).get("method", "otu")).lower()
        if any(str(r.get("platform", "illumina")).lower() in SHORT_P for r in amp):
            tools.append("dada2" if method == "asv" else "vsearch")
        if any(str(r.get("platform", "illumina")).lower() in (ONT_P | PB_P) for r in amp):
            tools.append("emu")
    tools.append("snakemake")
    # de-dupe, keep order
    seen, out = set(), []
    for t in tools:
        if t not in seen:
            seen.add(t)
            out.append(t)
    return out


def _scan_for_version(text: str) -> str:
    """First short, non-usage line carrying a dotted version, or 'unknown'."""
    junk = ("usage", "unrecognized", "error", "command not found", "must specify", "recommend")
    for line in (text or "").strip().splitlines():
        line = line.strip()
        if (line and len(line) < 160 and re.search(r"\d+\.\d+", line)
                and not any(j in line.lower() for j in junk)):
            return line
    return "unknown"


def tool_versions(tools: List[str]) -> Dict[str, str]:
    """Best-effort version capture; records 'not found' if a tool isn't on PATH.

    Prefers a registry's ``version_probe`` (the canonical command for that tool) when present —
    several bioinformatics tools print their version only via an idiosyncratic command (e.g.
    ``metabat2 --help`` carries ``version 2:2.18``), so the generic ``--version`` probe misreads
    or misses them. Falls back to trying ``--version``/``version``/``-v``/``-V``.
    """
    versions = {}
    for tool in tools:
        reg = registry.load_registry(tool) if tool in registry.list_tools() else {}
        cmd = reg.get("command", tool)
        exe = shutil.which(cmd)
        if not exe:
            versions[tool] = "not found on PATH"
            continue
        ver = "unknown"
        # 1) registry-declared canonical probe, if any
        probe = reg.get("version_probe")
        if probe:
            try:
                p = subprocess.run(shlex.split(probe), capture_output=True, text=True, timeout=20)
                ver = _scan_for_version((p.stdout or "") + "\n" + (p.stderr or ""))
            except (subprocess.SubprocessError, OSError):
                ver = "unknown"
        # 2) generic fallbacks
        if ver == "unknown":
            for flag in ("--version", "version", "-v", "-V"):
                try:
                    p = subprocess.run([cmd, flag], capture_output=True, text=True, timeout=15)
                    ver = _scan_for_version((p.stdout or "") + "\n" + (p.stderr or ""))
                    if ver != "unknown":
                        break
                except (subprocess.SubprocessError, OSError):
                    continue
        versions[tool] = ver
    return versions


def db_info(path: Optional[str]) -> Dict[str, Any]:
    if not path or not os.path.exists(path):
        return {"path": path, "present": False}
    info: Dict[str, Any] = {"path": os.path.abspath(path), "present": True}
    # hash the small kraken2 db hash header if available, else just record size
    hashfile = os.path.join(path, "hash.k2d")
    if os.path.isfile(hashfile):
        h = hashlib.md5()
        with open(hashfile, "rb") as fh:
            h.update(fh.read(1 << 20))  # first 1 MB is enough to fingerprint
        info["hash_k2d_md5_1mb"] = h.hexdigest()
    if os.path.isdir(path):
        total = sum(
            os.path.getsize(os.path.join(dp, f))
            for dp, _, fs in os.walk(path) for f in fs
        )
        info["size_bytes"] = total
    # If the DB was produced by db.build, carry its provenance manifest into the report so
    # the Methods/manifest record how (and from what) the reference DB was constructed.
    manifest = os.path.join(path, ".metagx_db.json") if os.path.isdir(path) else None
    if manifest and os.path.isfile(manifest):
        try:
            with open(manifest) as fh:
                info["build"] = json.load(fh)
        except (OSError, ValueError):
            pass
    return info


def _rendered_commands(cfg: Dict[str, Any]) -> Dict[str, List[str]]:
    """Re-render the exact command lines per tool (deterministic, registry-driven)."""
    cmds: Dict[str, List[str]] = {}
    mods = cfg.get("modules", {})
    if mods.get("classify"):
        base = dict(cfg.get("kraken2", {}))
        sweep = cfg.get("sweep")
        vals = sweep["values"] if sweep else [base.get("confidence", 0.0)]
        param = sweep["param"] if sweep else "confidence"
        sample = "<sample>"
        ex = []
        for v in vals:
            params = dict(base)
            params[param] = v
            args = registry.render_args(
                "kraken2", params,
                managed={"db": cfg["db"]["kraken2"], "threads": cfg.get("threads", 8),
                         "report": f"{sample}.{param}_{v}.kreport", "paired": True,
                         "gzip_compressed": True},
            )
            ex.append("kraken2 " + " ".join(args) + " <reads>")
        cmds["kraken2"] = ex
    if mods.get("abundance"):
        args = registry.render_args(
            "bracken", dict(cfg.get("bracken", {})),
            managed={"db": cfg["db"].get("bracken") or cfg["db"]["kraken2"],
                     "input": "<kreport>", "output": "<out>", "report_out": "<breport>"},
        )
        cmds["bracken"] = ["bracken " + " ".join(args)]
    if mods.get("qc"):
        args = registry.render_args(
            "fastp", dict(cfg.get("fastp", {})),
            managed={"thread": cfg.get("threads", 8), "in1": "<R1>", "out1": "<out_R1>"},
        )
        cmds["fastp"] = ["fastp " + " ".join(args)]
    if mods.get("assembly"):
        args = registry.render_args("megahit", dict(cfg.get("megahit", {})),
                                    managed={"threads": cfg.get("threads", 8)})
        cmds["megahit"] = ["megahit <reads> -o <out> " + " ".join(args)]
    if mods.get("binning"):
        args = registry.render_args("metabat2", dict(cfg.get("metabat2", {})),
                                    managed={"threads": cfg.get("threads", 8)})
        cmds["metabat2"] = ["metabat2 -i <contigs> -a <depth> -o <bins> " + " ".join(args)]
    return cmds


def qc_metrics(outdir: str) -> Dict[str, Any]:
    out = {}
    for jpath in sorted(glob.glob(os.path.join(outdir, "qc", "*.fastp.json"))):
        try:
            with open(jpath) as fh:
                j = json.load(fh)
            sample = os.path.basename(jpath)[: -len(".fastp.json")]
            out[sample] = {
                "reads_before": j.get("summary", {}).get("before_filtering", {}).get("total_reads"),
                "reads_after": j.get("summary", {}).get("after_filtering", {}).get("total_reads"),
            }
        except (json.JSONDecodeError, OSError):
            continue
    return out


def classification_metrics(outdir: str) -> Dict[str, Any]:
    """Percent classified per kraken2 report (100 - unclassified%)."""
    out = {}
    for rpath in sorted(glob.glob(os.path.join(outdir, "kraken2", "*.kreport"))):
        try:
            unclassified = 0.0
            with open(rpath) as fh:
                for line in fh:
                    cols = line.rstrip("\n").split("\t")
                    if len(cols) >= 4 and cols[-3].strip() == "U":
                        unclassified = float(cols[0])
                        break
            name = os.path.basename(rpath)[: -len(".kreport")]
            out[name] = {"percent_classified": round(100.0 - unclassified, 3)}
        except (ValueError, OSError):
            continue
    return out


# --------------------------------------------------------------------------- #
# Manifest + Methods                                                           #
# --------------------------------------------------------------------------- #
def build_manifest(cfg: Dict[str, Any]) -> Dict[str, Any]:
    outdir = _outdir(cfg)
    tools = active_tools(cfg)
    return {
        "metagx_version": __version__,
        "generated_at": _dt.datetime.now().astimezone().isoformat(timespec="seconds"),
        "project": cfg.get("project"),
        "preset": cfg.get("preset"),
        "config": cfg,
        "tool_versions": tool_versions(tools),
        "databases": {
            "kraken2": db_info(cfg.get("db", {}).get("kraken2")),
            "bracken": db_info(cfg.get("db", {}).get("bracken")),
        },
        "commands": _rendered_commands(cfg),
        "metrics": {
            "qc": qc_metrics(outdir),
            "classification": classification_metrics(outdir),
        },
    }


def methods_paragraph(cfg: Dict[str, Any], manifest: Dict[str, Any]) -> str:
    mods = cfg.get("modules", {})
    ver = manifest["tool_versions"]
    parts: List[str] = []

    def v(tool):
        val = ver.get(tool, "")
        if not val or val == "unknown" or "not found" in val:
            return ""  # only cite a real, captured version
        return f" (v{val})" if val[0].isdigit() else f" ({val})"

    plats = _platforms(cfg)
    amp = [r for r in _records(cfg) if str(r.get("library", "wgs")).lower() == "amplicon"]
    if mods.get("qc"):
        if plats & SHORT_P:
            fp = cfg.get("fastp", {})
            bits = []
            if "qualified_quality_phred" in fp:
                bits.append(f"a per-base quality threshold of Q{fp['qualified_quality_phred']}")
            if "length_required" in fp:
                bits.append(f"a minimum length of {fp['length_required']} bp")
            if fp.get("cut_right"):
                bits.append("3' sliding-window quality trimming")
            detail = (" with " + ", ".join(bits)) if bits else ""
            parts.append(
                f"Short reads were quality-filtered and adapter-trimmed with fastp{v('fastp')}{detail}."
            )
        if plats & ONT_P:
            parts.append(
                f"Nanopore reads had adapters removed with Porechop_ABI{v('porechop_abi')} "
                f"and were quality/length filtered with chopper{v('chopper')}."
            )
        if plats & PB_P:
            parts.append(f"PacBio reads were length/quality filtered with chopper{v('chopper')}.")

    if mods.get("classify"):
        kr = cfg.get("kraken2", {})
        sweep = cfg.get("sweep")
        if sweep:
            conf = ", ".join(str(x) for x in sweep["values"])
            cstr = f"across {sweep['param']} thresholds of {conf}"
        else:
            cstr = f"at a confidence threshold of {kr.get('confidence', 0.0)}"
        extra = []
        if "minimum_hit_groups" in kr:
            extra.append(f"a minimum of {kr['minimum_hit_groups']} hit groups")
        if kr.get("minimum_base_quality"):
            extra.append(f"a minimum base quality of {kr['minimum_base_quality']}")
        estr = (" requiring " + " and ".join(extra)) if extra else ""
        dbp = os.path.basename(str(cfg.get("db", {}).get("kraken2", "the reference database")))
        parts.append(
            f"Taxonomic classification was performed with Kraken 2{v('kraken2')} against "
            f"{dbp} {cstr}{estr}."
        )

    if mods.get("classify_consensus"):
        clf = str(cfg.get("consensus", {}).get("classifier", "metaphlan")).lower()
        if clf == "kaiju":
            parts.append(
                f"Classifications were cross-checked against Kaiju{v('kaiju')} (translated-"
                f"protein search); species-level concordance between the two methods is "
                f"reported per sample."
            )
        else:
            parts.append(
                f"Classifications were cross-checked against MetaPhlAn{v('metaphlan')} "
                f"(clade-specific markers); species-level concordance between the two "
                f"methods is reported per sample."
            )

    if mods.get("abundance"):
        br = cfg.get("bracken", {})
        parts.append(
            f"Species abundances were re-estimated with Bracken{v('bracken')} at the "
            f"{br.get('level', 'S')} level (read length {br.get('read_length', 100)}, "
            f"minimum {br.get('threshold', 0)} reads)."
        )

    if mods.get("assembly"):
        asm = str(cfg.get("assembly", {}).get("assembler", "megahit")).lower()
        if plats & SHORT_P:
            if asm == "metaspades":
                ms = cfg.get("metaspades", {})
                parts.append(
                    f"Short reads were assembled de novo with metaSPAdes{v('metaspades')} "
                    f"(SPAdes --meta, memory cap {ms.get('memory_gb', 250)} GB); where long "
                    f"reads accompanied a sample, a hybrid assembly was performed."
                )
            else:
                mh = cfg.get("megahit", {})
                pre = f" using the {mh['presets']} preset" if mh.get("presets") else ""
                parts.append(
                    f"Short reads were assembled de novo with MEGAHIT{v('megahit')}{pre} "
                    f"(minimum contig length {mh.get('min_contig_len', 200)} bp)."
                )
        if plats & (ONT_P | PB_P):
            fl = cfg.get("flye", {})
            mode = " in metagenome mode (metaFlye)" if fl.get("meta") else ""
            parts.append(f"Long reads were assembled de novo with Flye{v('flye')}{mode}.")
    if mods.get("binning"):
        mb = cfg.get("metabat2", {})
        parts.append(
            f"Reads were mapped back to contigs with minimap2{v('minimap2')} and "
            f"sorted with SAMtools{v('samtools')}; contigs were binned into MAGs with "
            f"MetaBAT 2{v('metabat2')} (minimum contig {mb.get('min_contig', 2500)} bp)."
        )
    if mods.get("bin_refinement"):
        dr = cfg.get("drep", {})
        parts.append(
            f"Bins were additionally recovered with MaxBin2{v('maxbin2')} and "
            f"CONCOCT{v('concoct')}; the three bin sets were reconciled into a consensus "
            f"per sample with DAS_Tool{v('das_tool')}, and the refined MAGs were "
            f"dereplicated across samples with dRep{v('drep')} (≥{dr.get('completeness', 75)}% "
            f"completeness, ≤{dr.get('contamination', 25)}% contamination, "
            f"{dr.get('s_ani', 0.95)} ANI)."
        )

    if mods.get("functional"):
        parts.append(
            f"Functional potential was profiled from reads with HUMAnN{v('humann')} "
            f"(MetaPhlAn{v('metaphlan')} taxonomic prescreen), yielding gene-family and "
            f"MetaCyc pathway abundances."
        )
        if mods.get("assembly"):
            parts.append(
                f"Antimicrobial-resistance and virulence determinants were screened on "
                f"assembled contigs with AMRFinderPlus{v('amrfinderplus')} and "
                f"ABRicate{v('abricate')}."
            )
        if mods.get("binning"):
            parts.append(
                f"Recovered bins were annotated with Bakta{v('bakta')} and assigned "
                f"orthology-based functions with eggNOG-mapper{v('eggnog')}."
            )

    if mods.get("aggregate"):
        parts.append(
            f"Per-sample QC and classification reports were aggregated with "
            f"MultiQC{v('multiqc')}, and community composition was visualized interactively "
            f"with Krona{v('krona')}."
        )
    if mods.get("strain"):
        parts.append(
            f"Within-population microdiversity (SNVs, nucleotide diversity) was profiled with "
            f"inStrain{v('instrain')} from the reads mapped to the assembly."
        )
    if mods.get("damage"):
        parts.append(
            f"Ancient-DNA authenticity was assessed with mapDamage2{v('mapdamage')}, "
            f"quantifying 5' C→T and 3' G→A post-mortem deamination at fragment ends."
        )
    if mods.get("decontam"):
        parts.append(
            "Reagent/laboratory contaminants were identified from the negative controls by a "
            "prevalence test and removed from the reported abundances."
        )
    if mods.get("differential"):
        d = cfg.get("differential", {})
        parts.append(
            f"Differential abundance between sample groups was tested on centred-log-ratio "
            f"transformed Bracken counts with a two-sided permutation test "
            f"({d.get('n_permutations', 999)} permutations), controlling the false-discovery "
            f"rate at {d.get('fdr', 0.05)} (Benjamini-Hochberg)."
        )
    if mods.get("bgc") and mods.get("assembly"):
        parts.append(
            f"Biosynthetic gene clusters encoding secondary metabolites were identified on "
            f"the assembled contigs with antiSMASH{v('antismash')}."
        )
    if amp:
        method = str(cfg.get("amplicon", {}).get("method", "otu")).lower()
        parts.append(
            f"Amplicon reads had primers removed with Cutadapt{v('cutadapt')}; short-read "
            + (f"libraries were denoised into amplicon sequence variants with DADA2{v('dada2')}"
               if method == "asv" else
               f"libraries were clustered into OTUs with VSEARCH{v('vsearch')}")
            + f", and long-read libraries were profiled with Emu{v('emu')}."
        )

    parts.append(f"The workflow was orchestrated with Snakemake{v('snakemake')}.")
    return " ".join(parts)


def citations_for(cfg: Dict[str, Any]) -> List[str]:
    return [CITATIONS[t] for t in active_tools(cfg) if t in CITATIONS]


# --------------------------------------------------------------------------- #
# Report rendering                                                             #
# --------------------------------------------------------------------------- #
def _param_table_md(cfg: Dict[str, Any]) -> str:
    rows = ["| Tool | Parameter | Value | Meaning |", "|---|---|---|---|"]
    for tool in registry.list_tools():
        section = cfg.get(tool)
        if not section:
            continue
        params = registry.load_registry(tool)["params"]
        for name, val in section.items():
            desc = " ".join(str(params.get(name, {}).get("question", "")).split())
            desc = (desc[:90] + "…") if len(desc) > 90 else desc
            rows.append(f"| {tool} | {name} | `{val}` | {desc} |")
    return "\n".join(rows) if len(rows) > 2 else "_No tool parameters customized._"


def _abundance_table_md(outdir: str, top_n: int = 15) -> str:
    path = os.path.join(outdir, "summary", "bracken_combined.tsv")
    if not os.path.isfile(path):
        return "_No Bracken abundance table found._"
    with open(path) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if not rows:
        return "_Bracken table is empty._"
    rows.sort(key=lambda r: float(r.get("fraction_total_reads", 0) or 0), reverse=True)
    head = ["sample", "label", "name", "new_est_reads", "fraction_total_reads"]
    head = [h for h in head if h in rows[0]]
    out = ["| " + " | ".join(head) + " |", "|" + "---|" * len(head)]
    for r in rows[:top_n]:
        out.append("| " + " | ".join(str(r.get(h, "")) for h in head) + " |")
    return "\n".join(out)


def _reconcile_md(outdir: str) -> str:
    """Summarize reconcile.json files: concordance counts + per-taxon read% vs contig%."""
    out = []
    for jpath in sorted(glob.glob(os.path.join(outdir, "reconcile", "*.reconcile.json"))):
        try:
            with open(jpath) as fh:
                s = json.load(fh)
        except (OSError, json.JSONDecodeError):
            continue
        c = s.get("taxa_concordance", {})
        out.append(
            f"**{s.get('sample')}** — {s.get('n_contigs_classified')}/{s.get('n_contigs')} "
            f"contigs classified, {s.get('assembled_bp')} bp assembled. Taxa concordance: "
            f"{c.get('both',0)} in both reads & contigs, {c.get('reads_only',0)} reads-only, "
            f"{c.get('contigs_only',0)} contigs-only; {s.get('n_flags',0)} discordance flag(s).\n"
        )
        cc = s.get("cat_cross_check")
        if cc:
            out.append(
                f"Contig taxonomy cross-check (kraken2 k-mer LCA vs CAT per-ORF voting): "
                f"{cc['kraken2_cat_agree']}/{cc['n_cat_classified']} agree, "
                f"{cc['kraken2_cat_conflict']} conflict.\n"
            )
        rows = s.get("top_taxa", [])[:10]
        if rows:
            out.append("| taxon | read % | contig cov-wt % | n contigs | concordance |\n|---|---|---|---|---|")
            for r in rows:
                out.append(f"| {r['taxon']} | {r['read_pct']} | {r['cov_weighted_pct']} | "
                           f"{r['n_contigs']} | {r['concordance']} |")
            out.append("")
    for jpath in sorted(glob.glob(os.path.join(outdir, "reconcile", "*.read_accuracy.json"))):
        try:
            with open(jpath) as fh:
                a = json.load(fh)
        except (OSError, json.JSONDecodeError):
            continue
        r = a.get("rates_on_classified_contigs", {})
        out.append(
            f"Read-level accuracy vs contig calls (**{a.get('sample')}**, "
            f"{a.get('reads_on_classified_contigs')} reads on classified contigs): "
            f"{r.get('concordant_pct')}% concordant, {r.get('discordant_pct')}% misclassified, "
            f"{r.get('read_unclassified_pct')}% unclassified.\n"
        )
    out.append("_read % = abundance from reads; contig cov-wt % = contig length×depth "
               "(taxonomy from the consensus). Divergence between them, or reads-only / "
               "contigs-only taxa, are flagged in `reconcile/*.flags.tsv`._\n")
    return "\n".join(out) if out else "_No reconciliation output found._\n"


def render_report_md(cfg: Dict[str, Any], manifest: Dict[str, Any]) -> str:
    outdir = _outdir(cfg)
    figs = sorted(glob.glob(os.path.join(outdir, "summary", "*.heatmap.png")))
    L: List[str] = []
    L.append(f"# Metagenomics report — {cfg.get('project', 'run')}\n")
    L.append(f"_Generated by metagx {manifest['metagx_version']} on {manifest['generated_at']}._\n")
    if cfg.get("preset"):
        L.append(f"**Preset:** `{cfg['preset']}`\n")

    L.append("## Methods\n")
    L.append(methods_paragraph(cfg, manifest) + "\n")

    L.append("## Software versions\n")
    L.append("| Tool | Version |\n|---|---|")
    for t, ver in manifest["tool_versions"].items():
        L.append(f"| {t} | {ver} |")
    L.append("")

    L.append("## Parameters\n")
    L.append(_param_table_md(cfg) + "\n")

    metrics = manifest.get("metrics", {})
    if metrics.get("classification"):
        L.append("## Classification metrics\n")
        L.append("| Sample/threshold | % classified |\n|---|---|")
        for k, m in metrics["classification"].items():
            L.append(f"| {k} | {m.get('percent_classified')} |")
        L.append("")

    if figs:
        L.append("## Figures\n")
        for f in figs:
            rel = os.path.relpath(f, outdir)
            L.append(f"![{os.path.basename(f)}]({rel})\n")

    if cfg.get("modules", {}).get("abundance"):
        L.append("## Top abundances (Bracken)\n")
        L.append(_abundance_table_md(outdir) + "\n")

    if cfg.get("modules", {}).get("reconcile"):
        L.append("## Read vs contig reconciliation\n")
        L.append(_reconcile_md(outdir))

    L.append("## References\n")
    for i, c in enumerate(citations_for(cfg), 1):
        L.append(f"{i}. {c}")
    L.append("")

    L.append("## Reproducibility\n")
    L.append("Exact command templates (managed I/O shown as placeholders):\n")
    L.append("```text")
    for tool, cmds in manifest.get("commands", {}).items():
        for c in cmds:
            L.append(c)
    L.append("```")
    return "\n".join(L)


def _md_to_latex(md_path: str, tex_path: str) -> bool:
    if not shutil.which("pandoc"):
        return False
    p = subprocess.run(["pandoc", md_path, "-o", tex_path], capture_output=True, text=True)
    return p.returncode == 0


def _md_to_pdf(md_path: str, pdf_path: str) -> bool:
    if not shutil.which("pandoc"):
        return False
    p = subprocess.run(["pandoc", md_path, "-o", pdf_path], capture_output=True, text=True)
    return p.returncode == 0


def generate(cfg: Dict[str, Any], fmt: str = "md") -> Dict[str, Any]:
    """Build manifest + Methods + report. ``fmt`` in {md, latex, pdf}.

    Returns a dict of written paths and the Methods text (for 'Copy as Methods').
    """
    outdir = _outdir(cfg)
    repdir = os.path.join(outdir, "report")
    os.makedirs(repdir, exist_ok=True)

    manifest = build_manifest(cfg)
    methods = methods_paragraph(cfg, manifest)
    report_md = render_report_md(cfg, manifest)

    written = {}
    mpath = os.path.join(repdir, "manifest.json")
    with open(mpath, "w") as fh:
        json.dump(manifest, fh, indent=2)
    written["manifest"] = mpath

    methpath = os.path.join(repdir, "methods.md")
    with open(methpath, "w") as fh:
        fh.write(methods + "\n\n## References\n\n")
        for i, c in enumerate(citations_for(cfg), 1):
            fh.write(f"{i}. {c}\n")
    written["methods"] = methpath

    rmd = os.path.join(repdir, "report.md")
    with open(rmd, "w") as fh:
        fh.write(report_md)
    written["report_md"] = rmd

    if fmt == "latex":
        tex = os.path.join(repdir, "report.tex")
        written["report_tex" if _md_to_latex(rmd, tex) else "latex_error"] = (
            tex if os.path.exists(tex) else "pandoc not available"
        )
    elif fmt == "pdf":
        pdf = os.path.join(repdir, "report.pdf")
        written["report_pdf" if _md_to_pdf(rmd, pdf) else "pdf_error"] = (
            pdf if os.path.exists(pdf) else "pandoc (+LaTeX) not available"
        )

    return {"paths": written, "methods": methods}
