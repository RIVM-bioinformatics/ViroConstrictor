"""Provide simple filesystem and JSON helpers for platform workflows."""

import json
from pathlib import Path
from typing import Any

import yaml

from container_manager.src.version import REPO_ROOT

MANIFESTS_DIR = REPO_ROOT / "container_manager" / "manifests"


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


def ensure_safe_path(
    path: Path,
    suffix: str,
    base_dir: Path = REPO_ROOT,
    must_exist: bool = False,
) -> Path:
    """Validate and resolve a path relative to a trusted base directory."""
    candidate = path.expanduser()
    resolved_path = candidate.resolve(strict=must_exist)
    resolved_base = base_dir.resolve(strict=True)

    try:
        resolved_path.relative_to(resolved_base)
    except ValueError as exc:
        raise ValueError(f"Unsafe path detected: {path} is outside {resolved_base}") from exc

    if resolved_path.suffix.lower() != suffix.lower():
        raise ValueError(f"Unsafe path detected: {path} does not end with {suffix}")
    return resolved_path


def ensure_safe_manifest_path(path: Path, must_exist: bool = False) -> Path:
    """Ensure a manifest path is a JSON file located under container_manager/manifests."""
    return ensure_safe_path(path=path, suffix=".json", base_dir=MANIFESTS_DIR, must_exist=must_exist)
