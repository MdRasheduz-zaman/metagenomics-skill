# Reproducible environments for metagx

Three levels of reproducibility, from convenient to bit-exact.

## 1. Per-rule conda envs (default, convenient)

The core tools come from `environment.yml`; every heavy optional tool has an isolated env
under `workflow/envs/`. Snakemake provisions them on first use:

```bash
metagx run --config config.yaml --use-conda
```

Versions are floor-pinned (`>=`), so a fresh solve can drift over time. Use a lock file or a
container when you need the *same* versions months later.

## 2. conda-lock (bit-level, no container)

[conda-lock](https://github.com/conda/conda-lock) resolves `environment.yml` (and the
`workflow/envs/*.yaml`) into a fully pinned, multi-platform lock file you commit to the repo:

```bash
pip install conda-lock

# core env -> conda-lock.yml (linux-64)
conda-lock lock -f environment.yml -p linux-64 --lockfile conda-lock.yml

# (optional) lock every per-rule env too
for e in workflow/envs/*.yaml; do
  conda-lock lock -f "$e" -p linux-64 --lockfile "containers/locks/$(basename "${e%.yaml}").conda-lock.yml"
done

# recreate the exact env later
conda-lock install --name metagx-core conda-lock.yml
```

Commit `conda-lock.yml` (and `containers/locks/*`) so a rebuild months later resolves to the
identical builds. Regenerate when you bump a tool in `environment.yml` or a `workflow/envs/*`.

## 3. Container (bit-exact, portable)

The repo `Dockerfile` builds an image with the metagx package + the core env:

```bash
docker build -t metagx:latest .
docker run --rm -v "$PWD:/data" metagx:latest \
    metagx run --config /data/config.yaml --use-conda
```

For HPC without Docker, convert to Apptainer/Singularity:

```bash
apptainer build metagx.sif docker-daemon://metagx:latest
apptainer exec --bind "$PWD:/data" metagx.sif metagx run --config /data/config.yaml
```

For the most reproducible image, point the `Dockerfile` at the committed `conda-lock.yml`
instead of `environment.yml` (`micromamba create -n base --file conda-lock.yml`).

> Note: per-rule Snakemake `container:` directives (one image per tool) are intentionally
> **deferred** — the per-rule conda envs + a single core image cover the common cases without
> maintaining a registry of per-tool images.
