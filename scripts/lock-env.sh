#!/usr/bin/env bash
# Regenerate the conda lockfile from environment.yml.
#
# environment.yml uses `>=` floors, which is right for "install me a working stack" but wrong
# for "reproduce the exact result I published": bioconda drifts, and two installs a month apart
# can resolve to different tool versions. conda-lock pins every transitive package to an exact
# build hash per platform, so a build from the lock is bit-for-bit reproducible. This is what
# the Dockerfile installs when conda-lock.yml is present.
#
# Run this after editing environment.yml, then commit the regenerated conda-lock.yml.
#
#   bash scripts/lock-env.sh
#
# Requires conda-lock (`pip install conda-lock` or `mamba install -c conda-forge conda-lock`).
set -euo pipefail
cd "$(dirname "$0")/.."

if ! command -v conda-lock >/dev/null 2>&1; then
  echo "conda-lock not found. Install it with:  pip install conda-lock" >&2
  exit 1
fi

# Lock linux-64 — the CI/Docker/HPC target, where bioconda has native builds for the whole
# stack so the env solves cleanly and reproducibly.
conda-lock lock \
  --file environment.yml \
  --platform linux-64 \
  --lockfile conda-lock.yml

echo "Wrote conda-lock.yml (linux-64) — commit it alongside environment.yml."
