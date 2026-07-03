#!/usr/bin/env bash
# ============================================================================
# Full test matrix for metagx. The heavy/real-execution tests are gated by
# ENVIRONMENT, not by code — because the tool envs deliberately COLLIDE:
# metagx-bio ships samtools >=1.18, while metagx-amr (abricate) hard-pins
# samtools 0.1.x and needs its own perl. You therefore CANNOT flatten every
# tool onto one PATH — the colliders run in their own env. That is the whole
# point of per-rule `--use-conda` isolation; this script mirrors it for tests.
#
# Each tier launches pytest with the .venv python EXPLICITLY (so metagx +
# snakemake always come from .venv) while putting the tool env on PATH for the
# subprocesses the tests shell out to.
#
#   bash scripts/run_all_tests.sh              # base + bio + amr + net
#   VIRAL=1 bash scripts/run_all_tests.sh      # also the ~0.6 GB viral fetch e2e
#
# Override env locations:  BIO=/path/to/metagx-bio/bin  AMR=/path/to/metagx-amr/bin
# ============================================================================
set -uo pipefail
REPO="$(cd "$(dirname "$0")/.." && pwd)"
PY="$REPO/.venv/bin/python"
BIO="${BIO:-$HOME/miniforge3/envs/metagx-bio/bin}"
AMR="${AMR:-$HOME/miniforge3/envs/metagx-amr/bin}"
VIRAL="${VIRAL:-0}"
rc=0

tier() { echo; echo "==================== $1 ===================="; }

# 1. BASE — pure-Python + workflow dry-run, no bio tools. The "runs anywhere" gate.
tier "base (no bio tools on PATH)"
env -u PATH PATH="$(dirname "$PY"):/usr/bin:/bin" "$PY" -m pytest "$REPO" -q || rc=1

# 2. BIO — real execution against metagx-bio (kraken2/bracken/fastp/BLAST validate/...).
tier "bio env (real e2e)"
PATH="$PATH:$BIO" "$PY" -m pytest "$REPO" -q || rc=1

# 3. AMR — colliders (abricate+blastn). metagx-amr FIRST so abricate's perl resolves; .venv
#    python still drives snakemake via sys.executable, so metagx is unaffected.
tier "amr env (colliders, isolated)"
PATH="$AMR:$PATH" "$PY" -m pytest \
    "$REPO/tests/test_functional_amr.py" \
    "$REPO/tests/test_pipeline_e2e.py::test_amr_screening_on_provided_genome" -q || rc=1

# 4. NET — live DB-URL liveness checks (opt-in; hits the network).
tier "network (METAGX_NET_TESTS=1)"
METAGX_NET_TESTS=1 "$PY" -m pytest "$REPO/tests/test_dbfetch.py" -q || rc=1

# 5. VIRAL — full RefSeq viral DB fetch (~0.6 GB) + real classification. Opt-in (slow).
if [ "$VIRAL" = "1" ]; then
  tier "viral fetch e2e (~0.6 GB download)"
  PATH="$PATH:$BIO" METAGX_E2E_VIRAL_FETCH=1 \
    "$PY" -m pytest "$REPO/tests/test_viral_fetch_e2e.py" -q || rc=1
fi

echo; echo "==================== RESULT ===================="
[ "$rc" -eq 0 ] && echo "ALL TIERS PASSED" || echo "SOME TIER FAILED (rc=$rc)"
exit "$rc"
