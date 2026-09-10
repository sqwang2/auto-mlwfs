#!/usr/bin/env python3
"""Shared machine-readable workflow output helpers."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any


def resolve_config(root: Path, configured: str | None, default_name: str) -> Path:
    """Resolve bundled defaults independently of the calculation directory."""

    if configured is None:
        return Path(__file__).resolve().parent.parent / "configuration" / default_name
    path = Path(configured).expanduser()
    return path if path.is_absolute() else root / path


def agent_output_dir(requested: Path) -> Path:
    """Choose a non-conflicting output path for an agent-prepared stage."""

    if not requested.exists():
        return requested
    alternative = requested.with_name(f"{requested.name}_agent")
    if alternative.exists():
        raise FileExistsError(
            f"Both {requested} and {alternative} already exist. "
            "Stop and ask the user which directory to use or preserve."
        )
    return alternative


def write_stage_summary(workdir: Path, summary: dict[str, Any]) -> Path:
    """Atomically write the stage result consumed by an agent."""

    path = workdir / "stage_summary.json"
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    temporary.replace(path)
    return path
