"""Workflow-specific helpers for container execution settings."""

import os
from pathlib import Path
from typing import Any


def construct_container_bind_args(samples_dict: dict[str, dict[str, Any]]) -> str:
    """Construct deduplicated bind mount arguments from existing sample file paths."""
    paths = [f"{Path(os.path.dirname(os.path.realpath(__file__))).parent}/"]
    for nested_dict in samples_dict.values():
        paths.extend(os.path.dirname(value) for value in nested_dict.values() if isinstance(value, str) and os.path.exists(value))
    return " ".join(f"--bind {path}" for path in set(paths))
