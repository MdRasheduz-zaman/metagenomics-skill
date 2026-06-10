"""Guard the tool-less paste-in prompt against drift.

prompts/INTERVIEW.md is the self-contained surface for clients with no CLI/MCP (Ollama,
plain chat). Unlike the registry-driven surfaces it is hand-written prose, so it silently
drifts when modules/presets are added. These tests fail if it falls behind the canonical
module list or the bundled presets/schedulers — turning a manual doc into a checked one.
"""
import os

from metagx import config_builder, presets, schedulers

_PROMPT = os.path.join(os.path.dirname(__file__), "..", "prompts", "INTERVIEW.md")


def _text():
    with open(_PROMPT) as fh:
        return fh.read()


def test_prompt_mentions_every_module():
    text = _text().lower()
    missing = [m for m in config_builder.DEFAULT_MODULES if m.lower() not in text]
    assert not missing, f"INTERVIEW.md is missing modules: {missing}"


def test_prompt_mentions_every_preset():
    text = _text()
    missing = [p["name"] for p in presets.describe_presets() if p["name"] not in text]
    assert not missing, f"INTERVIEW.md is missing presets: {missing}"


def test_prompt_mentions_hpc_executors():
    text = _text().lower()
    # the HPC backends should be offered to cluster users
    for name in ("slurm", "sge", "pbs", "lsf"):
        assert name in text, f"INTERVIEW.md should mention the {name} executor"
    assert "--executor" in text
    # sanity: the scheduler registry and the prompt agree slurm exists
    assert "slurm" in schedulers.list_schedulers()
