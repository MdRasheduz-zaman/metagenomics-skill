# metagx core environment as a container, for bit-level reproducibility.
#
#   Build:  docker build -t metagx:latest .
#   Run:    docker run --rm -v "$PWD:/data" metagx:latest metagx run --config /data/config.yaml
#
# The image carries the metagx package + the core tools from environment.yml (kraken2,
# bracken, fastp, megahit, minimap2, samtools, metabat2). The heavier optional tools
# (HUMAnN, GTDB-Tk, DAS_Tool, mapDamage2, inStrain, ...) stay in their per-rule conda envs
# under workflow/envs/; run with `metagx run --use-conda` (mount a conda pkg cache to reuse
# them) or build a fat image by adding `micromamba install -f workflow/envs/<env>.yaml`.
#
# Base image pinned by digest (not just the :1.5.8 tag, which can be re-pushed) for a
# reproducible build. This is the multi-arch manifest-list digest, so amd64/arm64 still
# resolve correctly. Refresh with: docker buildx imagetools inspect mambaorg/micromamba:1.5.8
FROM mambaorg/micromamba:1.5.8@sha256:475730daef12ff9c0733e70092aeeefdf4c373a584c952dac3f7bdb739601990

# Solve the environment into base. Prefer the committed conda-lock file (exact, hash-pinned,
# reproducible) when present; otherwise solve environment.yml. Regenerate the lock with
# `bash scripts/lock-env.sh` after editing environment.yml.
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
COPY --chown=$MAMBA_USER:$MAMBA_USER conda-lock.ym[l] /tmp/
RUN if [ -f /tmp/conda-lock.yml ]; then \
        micromamba install -y -n base -f /tmp/conda-lock.yml ; \
    else \
        micromamba install -y -n base -f /tmp/environment.yml ; \
    fi && \
    micromamba clean --all --yes

# Install the metagx package itself.
ARG MAMBA_DOCKERFILE_ACTIVATE=1
COPY --chown=$MAMBA_USER:$MAMBA_USER . /opt/metagx
RUN pip install --no-cache-dir -e /opt/metagx

WORKDIR /data
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh"]
CMD ["metagx", "--help"]
