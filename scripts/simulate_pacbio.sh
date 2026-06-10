#!/usr/bin/env bash
# Simulate PacBio CLR and HiFi long reads from a reference FASTA for metagx smoke tests.
#
# Tools (TEST-ONLY — remove from metagx-bio when done):
#   pbsim3  — CLR via PBSIM quality model; multipass BAM for HiFi+ccs on Linux
#   pbccs   — CCS consensus from multipass BAM (Linux binary; skipped on macOS)
#   simlord — PacBio-like reads; used for HiFi on macOS and as CLR alternate
#
#   conda install -c bioconda pbsim3 pbccs simlord
#   conda remove -n metagx-bio pbsim3 pbccs simlord
#
# Usage:
#   bash scripts/simulate_pacbio.sh
#   bash scripts/simulate_pacbio.sh --depth 15 --kept 1500
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
REF="${REF:-$ROOT/data/genomes.fasta}"
OUT="${OUT:-$ROOT/data/pacbio_sim}"
DEPTH="${DEPTH:-5}"
SEED="${SEED:-42}"
KEPT="${KEPT:-2000}"
MODEL="${MODEL:-${CONDA_PREFIX:-}/data/QSHMM-RSII.model}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --depth) DEPTH="$2"; shift 2 ;;
    --kept) KEPT="$2"; shift 2 ;;
    --ref) REF="$2"; shift 2 ;;
    --out) OUT="$2"; shift 2 ;;
    -h|--help)
      sed -n '1,22p' "$0"
      exit 0
      ;;
    *) echo "Unknown arg: $1" >&2; exit 2 ;;
  esac
done

for cmd in pbsim seqtk gzip gunzip; do
  command -v "$cmd" >/dev/null || { echo "$cmd not on PATH — activate metagx-bio" >&2; exit 1; }
done
command -v simlord >/dev/null || { echo "simlord not on PATH — conda install -c bioconda simlord" >&2; exit 1; }

if [[ ! -f "$MODEL" ]]; then
  for candidate in \
    "$HOME/miniconda3/envs/metagx-bio/data/QSHMM-RSII.model" \
    "$(dirname "$(command -v pbsim)")/../data/QSHMM-RSII.model"; do
    if [[ -f "$candidate" ]]; then MODEL="$candidate"; break; fi
  done
fi
[[ -f "$MODEL" ]] || { echo "PBSIM QSHMM model not found (set MODEL=)" >&2; exit 1; }

ccs_runnable() {
  command -v ccs >/dev/null || return 1
  ccs --version >/dev/null 2>&1
}

mkdir -p "$OUT"
TMP_CLR="$OUT/_pbsim_clr_work"
TMP_HIFI="$OUT/_pbsim_hifi_work"
SPLIT_DIR="$OUT/_split_refs"
rm -rf "$TMP_CLR" "$TMP_HIFI" "$SPLIT_DIR"
mkdir -p "$TMP_CLR" "$TMP_HIFI" "$SPLIT_DIR"

# PBSIM3 can exit early on multi-record FASTA on some hosts; simulate one ref at a time.
python3 - "$REF" "$SPLIT_DIR" <<'PY'
import sys
from pathlib import Path

ref, out = Path(sys.argv[1]), Path(sys.argv[2])
acc, lines, idx = None, [], 0
for line in ref.open():
    if line.startswith(">"):
        if acc is not None:
            (out / f"ref_{idx:04d}.fasta").write_text("".join(lines))
            idx += 1
        acc = line[1:].split()[0]
        lines = [line]
    elif acc is not None:
        lines.append(line)
if acc is not None:
    (out / f"ref_{idx:04d}.fasta").write_text("".join(lines))
print(idx + 1, "references split")
PY

run_pbsim_per_ref() {
  local work="$1" prefix="$2" passes="$3"
  local i=0
  for one in "$SPLIT_DIR"/ref_*.fasta; do
    i=$((i + 1))
    (
      cd "$work"
      pbsim --strategy wgs --method qshmm --qshmm "$MODEL" \
        --depth "$DEPTH" --genome "$one" --prefix "${prefix}_${i}" --seed "$SEED" \
        --pass-num "$passes" --length-min 500 --length-max 25000 \
        > "${prefix}_${i}.log" 2>&1
    )
  done
}

echo "==> PBSIM3 CLR (pass-num 1, depth=$DEPTH) from $(basename "$REF")"
run_pbsim_per_ref "$TMP_CLR" clr 1

CLR_MERGED="$OUT/_clr_merged.fastq.gz"
echo "==> Merge PBSIM CLR shards"
shopt -s nullglob
CLR_SHARDS=( "$TMP_CLR"/clr_*.fq.gz )
if ((${#CLR_SHARDS[@]} == 0)); then
  echo "No PBSIM CLR output in $TMP_CLR (see ${TMP_CLR}/clr_*.log)" >&2
  exit 1
fi
gunzip -c "${CLR_SHARDS[@]}" | gzip -c > "$CLR_MERGED"
echo "==> Subsample CLR -> pacbio_clr.fastq.gz ($KEPT reads)"
seqtk sample -s"$SEED" "$CLR_MERGED" "$KEPT" | gzip -c > "$OUT/pacbio_clr.fastq.gz"
rm -f "$CLR_MERGED"
rm -rf "$TMP_CLR"

HIFI_OUT="$OUT/pacbio_hifi.fastq.gz"
if ccs_runnable; then
  echo "==> PBSIM3 multipass + CCS HiFi (depth=$DEPTH, pass-num=10)"
  run_pbsim_per_ref "$TMP_HIFI" hifi 10
  HIFI_PARTS=()
  for bam in "$TMP_HIFI"/hifi_*_*.bam; do
    base="$(basename "$bam" .bam)"
    ccs "$bam" "${base}.ccs.fastq.gz" --num-threads 2
    HIFI_PARTS+=( "${TMP_HIFI}/${base}.ccs.fastq.gz" )
  done
  gunzip -c "${HIFI_PARTS[@]}" | gzip -c > "$OUT/_hifi_merged.fastq.gz"
  seqtk sample -s"$SEED" "$OUT/_hifi_merged.fastq.gz" "$KEPT" | gzip -c > "$HIFI_OUT"
  rm -f "$OUT/_hifi_merged.fastq.gz"
  rm -rf "$TMP_HIFI"
  echo "    HiFi source: PBSIM3 + CCS"
else
  echo "==> CCS not runnable on this host — SimLoRD HiFi (coverage=$DEPTH, max-passes=10)"
  simlord --read-reference "$REF" --coverage "$DEPTH" --max-passes 10 \
    --gzip --no-sam "$OUT/_simlord_hifi" >/dev/null
  seqtk sample -s"$SEED" "$OUT/_simlord_hifi.fastq.gz" "$KEPT" | gzip -c > "$HIFI_OUT"
  rm -f "$OUT/_simlord_hifi.fastq.gz"
  echo "    HiFi source: SimLoRD multipass (use Linux + pbccs for PBSIM3+CCS)"
fi

echo "==> SimLoRD CLR alternate (coverage=$DEPTH, max-passes=1)"
simlord --read-reference "$REF" --coverage "$DEPTH" --max-passes 1 \
  --gzip --no-sam "$OUT/_simlord_clr" >/dev/null
seqtk sample -s"$SEED" "$OUT/_simlord_clr.fastq.gz" "$KEPT" | gzip -c > "$OUT/pacbio_clr_simlord.fastq.gz"
rm -f "$OUT/_simlord_clr.fastq.gz"
rm -rf "$SPLIT_DIR"

echo "==> Done. Outputs in $OUT"
ls -lh "$OUT"/*.fastq.gz
