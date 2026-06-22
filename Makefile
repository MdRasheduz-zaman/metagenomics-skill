# metagx developer targets. `make help` lists them.
.PHONY: help test dryrun e2e repro-ci lock

help:
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) \
		| sort | awk 'BEGIN{FS=":.*?## "}{printf "  %-12s %s\n", $$1, $$2}'

test: ## Full unit suite + workflow dry-run gate (no bio tools needed)
	pytest -q
	pytest -q tests/test_workflow_dryrun.py

dryrun: ## Workflow DAG / render_args gate only
	pytest -q tests/test_workflow_dryrun.py

e2e: ## Real kraken2+Bracken end-to-end (needs metagx-bio tools on PATH)
	pytest -q tests/test_pipeline_e2e.py tests/test_adna_e2e.py

# Reproduce the low-core CI environment that hid the thread>core DB-build bugs for four
# rounds. Forces a from-scratch DB build (METAGX_FORCE_DB_BUILD) under 2 visible CPUs so
# bracken-build/kmer2read_distr and kraken2-build hit the same constraint as the 2-core
# GitHub runner. On Linux this pins cores with taskset (real reproduction: without the
# dbbuild thread-clamp, kmer2read_distr would hard-fail here). Elsewhere it falls back to
# OMP_NUM_THREADS=2, which only constrains OpenMP — weaker, but still exercises the path.
repro-ci: ## Reproduce the 2-core CI DB-build path locally (catches thread>core regressions)
	@if command -v taskset >/dev/null 2>&1; then \
		echo ">> taskset -c 0,1 (real 2-core pinning)"; \
		OMP_NUM_THREADS=2 METAGX_FORCE_DB_BUILD=1 taskset -c 0,1 \
			pytest -q tests/test_pipeline_e2e.py::test_platform_classifies_correctly; \
	else \
		echo ">> taskset unavailable (non-Linux): OMP_NUM_THREADS=2 only — weaker reproduction"; \
		OMP_NUM_THREADS=2 METAGX_FORCE_DB_BUILD=1 \
			pytest -q tests/test_pipeline_e2e.py::test_platform_classifies_correctly; \
	fi

lock: ## Regenerate the conda lockfile
	bash scripts/lock-env.sh
