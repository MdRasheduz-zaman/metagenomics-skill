"""`metagx refresh`: propose registry updates from the installed binary.

The core (`plan_refresh`) is pure over an injected capture, so it runs with no tool installed —
same style as tests/test_toollock.py. A new binary flag must surface as a `_status: proposed`
draft, a registry flag absent from the binary must be reported as drift, and a version delta must
propose a tested_version bump.
"""
from metagx import refresh


def _cap(version, help_text, flags):
    """A fake sync_help.capture_help result (ok=True)."""
    return {"ok": True, "version": version, "help_text": help_text,
            "flags": [{"flag": f} for f in flags]}


_HELP = """\
Usage: kraken2 [options] <files>
  --confidence FLOAT      confidence score threshold
  --brand-new-flag INT    a brand new option not yet in the registry
"""


def test_new_binary_flag_becomes_a_proposed_stub():
    cap = _cap("Kraken version 2.17.1", _HELP, ["--confidence", "--brand-new-flag"])
    prop = refresh.plan_refresh("kraken2", capture=cap)
    assert prop["capture_ok"]
    assert "brand_new_flag" in prop["new_params"]
    stub = prop["new_params"]["brand_new_flag"]
    assert stub["_status"] == "proposed"        # inert until curated
    assert stub["flag"] == "--brand-new-flag"
    assert stub["type"] == "int"                # guessed from the INT metavar
    assert stub["question"] == "<LLM draft — REVIEW>"


def test_registry_flag_absent_from_binary_is_reported_as_drift():
    # help lists only --confidence, so every other curated kraken2 flag is "removed".
    cap = _cap("Kraken version 2.17.1", _HELP, ["--confidence", "--brand-new-flag"])
    prop = refresh.plan_refresh("kraken2", capture=cap)
    assert "--use-mpa-style" in prop["removed_flags"]


def test_version_delta_flags_a_bump():
    cap = _cap("Kraken version 2.18.0", _HELP, ["--confidence"])
    prop = refresh.plan_refresh("kraken2", capture=cap)
    assert prop["version"]["differs"] is True
    assert prop["version"]["tested"] == "2.17.1"
    assert "version drift" in prop["summary"]


def test_matching_version_no_delta():
    cap = _cap("Kraken version 2.17.1", _HELP, ["--confidence"])
    prop = refresh.plan_refresh("kraken2", capture=cap)
    assert prop["version"]["differs"] is False


def test_capture_failure_proposes_nothing():
    prop = refresh.plan_refresh("kraken2", capture={"ok": False, "error": "not on PATH"})
    assert prop["capture_ok"] is False
    assert prop["new_params"] == {} and prop["removed_flags"] == []


def test_write_proposal_uses_injected_writer_and_never_touches_registry():
    cap = _cap("Kraken version 2.18.0", _HELP, ["--confidence", "--brand-new-flag"])
    prop = refresh.plan_refresh("kraken2", capture=cap)
    written = {}
    path = refresh.write_proposal(prop, root="/somewhere",
                                  writer=lambda p, c: written.__setitem__(p, c))
    assert path.endswith("refresh/kraken2.proposed.yaml")
    assert "_status: proposed" in written[path]
    assert "brand_new_flag" in written[path]


def test_write_proposal_noop_when_nothing_new():
    cap = _cap("Kraken version 2.17.1", _HELP, ["--confidence"])  # no new flags
    prop = refresh.plan_refresh("kraken2", capture=cap)
    assert refresh.write_proposal(prop, root="/x", writer=lambda p, c: None) is None


def test_version_token_extraction():
    assert refresh.version_token("Kraken version 2.17.1") == "2.17.1"
    assert refresh.version_token("v1.0") == "1.0"
    assert refresh.version_token(None) is None
    assert refresh.version_token("no digits here") is None
