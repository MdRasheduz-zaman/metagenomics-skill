#!/usr/bin/env python3
"""Thin wrapper over `metagx.compare` — see `metagx compare --help`.

    PATH="$HOME/miniconda3/envs/metagx-bio/bin:$PATH" \
        .venv/bin/python scripts/compare_platforms.py [manifest.tsv] [outdir]

Equivalent to `metagx compare [--manifest ...] [--outdir ...]`. Kept so the experiment-08
reproducibility recipe (and any agent that shells scripts) keeps working.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from metagx.compare import main  # noqa: E402

if __name__ == "__main__":
    raise SystemExit(main())
