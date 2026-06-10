#!/usr/bin/env bash
# Simulate Illumina paired-end reads from a reference FASTA with wgsim, then subsample
# into study-sized FASTQ files for metagx smoke tests.
#
# TEST-ONLY TOOL: wgsim is for local rule/pipeline validation. Real projects use
# experimental FASTQ; remove wgsim from the metagx-bio env when disk is tight:
#   conda remove -n metagx-bio wgsim
#
# Requires: wgsim, seqtk, gzip (metagx-bio conda env recommended).
# Usage:
#   bash scripts/simulate_illumina_wgsim.sh
#   bash scripts/simulate_illumina_wgsim.sh --pairs 50000 --pe-kept 10000 --diff-kept 3000
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
REF="${REF:-$ROOT/data/genomes.fasta}"
OUT="${OUT:-$ROOT/data/illumina_sim}"
PAIRS="${PAIRS:-50000}"
READ_LEN="${READ_LEN:-150}"
INSERT="${INSERT:-350}"
INSERT_SD="${INSERT_SD:-50}"
ERROR_RATE="${ERROR_RATE:-0.02}"
SEED="${SEED:-42}"
PE_KEPT="${PE_KEPT:-10000}"
DIFF_KEPT="${DIFF_KEPT:-3000}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --pairs) PAIRS="$2"; shift 2 ;;
    --pe-kept) PE_KEPT="$2"; shift 2 ;;
    --diff-kept) DIFF_KEPT="$2"; shift 2 ;;
    --ref) REF="$2"; shift 2 ;;
    --out) OUT="$2"; shift 2 ;;
    -h|--help)
      sed -n '1,20p' "$0"
      exit 0
      ;;
    *) echo "Unknown arg: $1" >&2; exit 2 ;;
  esac
done

command -v wgsim >/dev/null || { echo "wgsim not on PATH — activate metagx-bio or: conda install -c bioconda wgsim" >&2; exit 1; }
command -v seqtk >/dev/null || { echo "seqtk not on PATH" >&2; exit 1; }
command -v python3 >/dev/null || { echo "python3 not on PATH" >&2; exit 1; }

# wgsim writes placeholder Q2 scores ('2'); kraken2 --minimum-base-quality needs real Illumina Q.
fix_illumina_quals() {
  python3 - "$1" "$2" <<'PY'
import sys
qchar = chr(33 + 30)  # Phred+33 Q30
path_in, path_out = sys.argv[1], sys.argv[2]
with open(path_in) as fin, open(path_out, "w") as fout:
    for i, line in enumerate(fin, 1):
        if i % 4 == 0:
            fout.write(qchar * len(line.rstrip("\n")) + "\n")
        else:
            fout.write(line)
PY
}

mkdir -p "$OUT"
RAW_R1="$OUT/_wgsim_R1.fastq"
RAW_R2="$OUT/_wgsim_R2.fastq"
FIX_R1="$OUT/_wgsim_R1.q30.fastq"
FIX_R2="$OUT/_wgsim_R2.q30.fastq"

echo "==> wgsim: $PAIRS PE pairs, ${READ_LEN}bp, from $(basename "$REF")"
wgsim -N "$PAIRS" -1 "$READ_LEN" -2 "$READ_LEN" \
  -d "$INSERT" -s "$INSERT_SD" -e "$ERROR_RATE" \
  -r 0.001 -R 0 -S "$SEED" \
  "$REF" "$RAW_R1" "$RAW_R2"

echo "==> Upgrade wgsim Q2 placeholders to Phred+33 Q30 (Illumina-like)"
fix_illumina_quals "$RAW_R1" "$FIX_R1"
fix_illumina_quals "$RAW_R2" "$FIX_R2"
rm -f "$RAW_R1" "$RAW_R2"
RAW_R1="$FIX_R1"
RAW_R2="$FIX_R2"

echo "==> PE smoke sample: seqtk sample $PE_KEPT pairs -> illumina_sim_R{1,2}.fastq.gz"
seqtk sample -s"$SEED" "$RAW_R1" "$PE_KEPT" | gzip -c > "$OUT/illumina_sim_R1.fastq.gz"
seqtk sample -s"$SEED" "$RAW_R2" "$PE_KEPT" | gzip -c > "$OUT/illumina_sim_R2.fastq.gz"

echo "==> Differential pseudo-samples (SE, R1 only, different seeds)"
for label in case1 case2 ctrl1 ctrl2; do
  s=$((SEED + $(echo -n "$label" | cksum | awk '{print $1 % 10000}')))
  seqtk sample -s"$s" "$RAW_R1" "$DIFF_KEPT" | gzip -c > "$OUT/${label}_R1.fastq.gz"
done

rm -f "$RAW_R1" "$RAW_R2"
echo "==> Done. Outputs in $OUT (qualities set to Q30 for kraken2 --minimum-base-quality)"
ls -lh "$OUT"/*.fastq.gz
