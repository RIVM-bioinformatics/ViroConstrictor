"""Provide simple filesystem and JSON helpers for platform workflows."""

import json
from pathlib import Path
from typing import Any

import yaml


def ensure_dir(path: Path) -> None:
    """Create a directory path recursively if it does not already exist."""
    path.mkdir(parents=True, exist_ok=True)


def write_json(path: Path, payload: dict[str, Any]) -> None:
    """Write a dictionary to a JSON file with stable indentation and newline."""
    ensure_dir(path.parent)
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def read_json(path: Path) -> dict[str, Any]:
    """Read and parse a JSON file into a dictionary."""
    return json.loads(path.read_text(encoding="utf-8"))


def load_yaml(path: Path) -> dict[str, Any]:
    """Load a YAML file and return its top-level mapping as a dictionary."""
    return yaml.safe_load(path.read_text(encoding="utf-8"))
