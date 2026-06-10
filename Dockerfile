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
# Pin to a specific micromamba/base tag + a conda-lock file (see containers/README.md) for a
# fully reproducible build.
FROM mambaorg/micromamba:1.5.8

# Solve the core environment into base.
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean --all --yes

# Install the metagx package itself.
ARG MAMBA_DOCKERFILE_ACTIVATE=1
COPY --chown=$MAMBA_USER:$MAMBA_USER . /opt/metagx
RUN pip install --no-cache-dir -e /opt/metagx

WORKDIR /data
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh"]
CMD ["metagx", "--help"]
