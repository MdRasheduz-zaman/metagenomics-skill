#!/usr/bin/env bash
# Assemble each sequencing platform's reads from the shared 30-genome reference into a
# common location so they can be compared side-by-side (scripts/compare_platforms.py).
#
#   short reads (Illumina) -> MEGAHIT      long reads (ONT/PacBio) -> Flye
#
# This mirrors metagx's platform-routed assembly (workflow/rules/assembly.smk). On real
# Linux/Docker hosts, prefer `metagx run --config <cfg>` with `modules.assembly: true`.
# This standalone driver exists because the bioconda MEGAHIT is an x86_64 binary that is
# unstable under Rosetta on arm64 macOS; here it is invoked with the Rosetta-safe flags
# (--no-hw-accel, constrained k-list, single thread) found to complete on this host.
#
# Usage: bash scripts/run_platform_assemblies.sh
set -uo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
export PATH="${HOME}/miniconda3/envs/metagx-bio/bin:${PATH}"
OUT="$ROOT/results/experiments/cross_platform_assembly"
THREADS="${THREADS:-4}"
mkdir -p "$OUT"

contig_count() { grep -c '^>' "$1" 2>/dev/null || echo 0; }

# --- Illumina: MEGAHIT (short-read), Rosetta-safe flags -----------------------------
ilmn_out="$OUT/illumina"
mkdir -p "$ilmn_out"
echo "==> Illumina / MEGAHIT (short-read)"
rm -rf "$ilmn_out/_mh"
if megahit --no-hw-accel --k-list 21,29,39 -t 1 \
     -1 "$ROOT/data/illumina_sim/illumina_sim_R1.fastq.gz" \
     -2 "$ROOT/data/illumina_sim/illumina_sim_R2.fastq.gz" \
     -o "$ilmn_out/_mh" --min-contig-len 200 > "$ilmn_out/megahit.log" 2>&1; then
  cp "$ilmn_out/_mh/final.contigs.fa" "$ilmn_out/contigs.fa"
  rm -rf "$ilmn_out/_mh"
  echo "    OK: $(contig_count "$ilmn_out/contigs.fa") contigs"
else
  echo "    FAILED (see $ilmn_out/megahit.log)"; : > "$ilmn_out/contigs.fa"
fi

# --- ONT: reuse the Flye assembly already produced by the metagx pipeline ------------
ont_out="$OUT/ont"
mkdir -p "$ont_out"
echo "==> ONT / Flye (long-read) — reuse metagx output"
ont_src="$ROOT/results/ont_sim/assembly/ont_sim/final.contigs.fa"
if [[ -s "$ont_src" ]]; then
  cp "$ont_src" "$ont_out/contigs.fa"
  echo "    OK: $(contig_count "$ont_out/contigs.fa") contigs (from $ont_src)"
else
  echo "    ONT contigs not found — run experiment 03 (ont assembly) first"; : > "$ont_out/contigs.fa"
fi

# --- PacBio HiFi + CLR: Flye (long-read, metagenome) ---------------------------------
run_flye() {
  local label="$1" flag="$2" reads="$3"
  local d="$OUT/$label"; mkdir -p "$d"
  echo "==> ${label} / Flye ${flag} (long-read)"
  rm -rf "$d/_flye"
  if flye "$flag" "$reads" --out-dir "$d/_flye" -t "$THREADS" --meta > "$d/flye.log" 2>&1; then
    cp "$d/_flye/assembly.fasta" "$d/contigs.fa"; rm -rf "$d/_flye"
    echo "    OK: $(contig_count "$d/contigs.fa") contigs"
  else
    # Flye aborts on too-low coverage ("No disjointigs were assembled"); record empty.
    echo "    NO ASSEMBLY: $(grep -m1 'ERROR' "$d/flye.log" | sed 's/.*ERROR: //')"
    : > "$d/contigs.fa"; rm -rf "$d/_flye"
  fi
}
run_flye pacbio_hifi --pacbio-hifi "$ROOT/data/pacbio_sim/pacbio_hifi.fastq.gz"
run_flye pacbio_clr  --pacbio-raw  "$ROOT/data/pacbio_sim/pacbio_clr_simlord.fastq.gz"

echo "==> Done. Contigs under $OUT/<platform>/contigs.fa"