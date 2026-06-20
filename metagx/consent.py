"""Tiny persisted consent store for capabilities that touch user data.

The only scope today is ``probe`` (local read profiling). Consent is remembered in
``.metagx/consent.json`` so we never re-nag, and the trust boundary is explicit: a stored
``local`` consent never implies permission to send anything off the machine — any future
off-box capability must ask again under its own scope.
"""

from __future__ import annotations

import json
import os
import time
from typing import Optional

STORE = os.path.join(".metagx", "consent.json")
VALID = {"local", "off"}


def _load() -> dict:
    try:
        with open(STORE) as fh:
            return json.load(fh)
    except (OSError, ValueError):
        return {}


def get(scope: str) -> Optional[str]:
    """Stored consent for ``scope`` (e.g. 'probe'), or None if never decided."""
    val = _load().get(scope, {}).get("value")
    return val if val in VALID else None


def set(scope: str, value: str) -> str:
    """Persist ``value`` ('local' | 'off') for ``scope``. Returns the value."""
    if value not in VALID:
        raise ValueError(f"consent value must be one of {sorted(VALID)}, got {value!r}")
    data = _load()
    data[scope] = {"value": value, "ts": int(time.time())}
    os.makedirs(os.path.dirname(STORE), exist_ok=True)
    with open(STORE, "w") as fh:
        json.dump(data, fh, indent=2)
    return value
