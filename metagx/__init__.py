"""metagx: schema-driven metagenomics pipeline core.

The parameter registries under ``metagx/parameters/*.yaml`` are the single source of
truth. Everything else (the interview, config validation, MCP tool schemas, the CLI,
and command-line construction inside the Snakefile) is derived from them.
"""

from .registry import (
    list_tools,
    load_registry,
    interview_spec,
    render_args,
    validate,
    ValidationError,
)

__all__ = [
    "list_tools",
    "load_registry",
    "interview_spec",
    "render_args",
    "validate",
    "ValidationError",
]

__version__ = "0.1.0"
