#!/usr/bin/env bash
# Install the bioinformatics tools on Apple Silicon (arm64) macOS, where bioconda has no
# native kraken2/bracken builds. We use an osx-64 (Rosetta) conda env for kraken2/fastp/
# seqtk, and build Bracken from source (its only osx-64 conda build is a broken placeholder).
#
# Verified on macOS 14 (Darwin 24) with miniconda + Homebrew. Re-runnable.
set -euo pipefail

ENV_NAME="${1:-metagx-bio}"
CONDA="${CONDA_EXE:-$HOME/miniconda3/bin/conda}"

echo "==> 1/4 osx-64 conda env '$ENV_NAME' (kraken2, fastp, seqtk)"
CONDA_SUBDIR=osx-64 "$CONDA" create -y -n "$ENV_NAME" -c conda-forge -c bioconda \
    python=3.11 kraken2 fastp seqtk

ENV_DIR="$("$CONDA" env list | awk -v n="$ENV_NAME" '$1==n{print $NF}')"
[ -n "$ENV_DIR" ] || { echo "could not locate env dir"; exit 1; }

echo "==> 2/4 libomp (for compiling Bracken's OpenMP helper with Apple clang)"
command -v brew >/dev/null || { echo "Homebrew required: https://brew.sh"; exit 1; }
brew list libomp >/dev/null 2>&1 || brew install libomp
LIBOMP="$(brew --prefix libomp)"

echo "==> 3/4 build Bracken from source"
SRC="$ENV_DIR/opt/Bracken"
rm -rf "$SRC"
git clone --depth 1 https://github.com/jenniferlu717/Bracken.git "$SRC"
make -C "$SRC/src" clean >/dev/null 2>&1 || true
make -C "$SRC/src" CXX=clang++ \
    CXXFLAGS="-c -g -O3 -std=c++11 -Xpreprocessor -fopenmp -I$LIBOMP/include" \
    LDFLAGS="-L$LIBOMP/lib -lomp -Wl,-rpath,$LIBOMP/lib"

# Upstream bug: in the branch that takes -w/--out-report, the est_abundance.py call is
# prefixed with a stray `echo`, so Bracken prints "complete" but writes nothing. Remove it.
perl -0pi -e 's/(\n\s*)echo (python \$DIR\/src\/est_abundance\.py -i \$\{INPUT\} \\\n\s*-o \$\{OUTPUT\} \\\n\s*--out-report)/$1$2/' "$SRC/bracken"

ln -sf "$SRC/bracken" "$ENV_DIR/bin/bracken"
ln -sf "$SRC/bracken-build" "$ENV_DIR/bin/bracken-build"

echo "==> 4/4 verify"
export PATH="$ENV_DIR/bin:$PATH"
for t in kraken2 kraken2-build bracken bracken-build fastp seqtk; do
    printf "  %-14s %s\n" "$t" "$(command -v "$t" || echo MISSING)"
done
echo
echo "Done. Put the env on PATH before running metagx, e.g.:"
echo "    export PATH=\"$ENV_DIR/bin:\$PATH\""
