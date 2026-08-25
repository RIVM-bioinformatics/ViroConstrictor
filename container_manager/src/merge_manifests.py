"""Merge container build manifests produced by matrix jobs into one manifest."""

import glob
import json
from pathlib import Path
from typing import Any

from container_manager.src.io import read_json


def merge_manifests(input_glob: str, output_path: Path) -> int:
    """Merge all manifest files matching input_glob into output_path."""
    manifest_paths = [Path(p) for p in sorted(glob.glob(input_glob))]
    if not manifest_paths:
        raise SystemExit(f"No manifest files matched pattern: {input_glob}")

    merged: dict[str, Any] | None = None
    items: list[dict[str, Any]] = []

    for manifest_path in manifest_paths:
        data = read_json(manifest_path)
        if merged is None:
            merged = {key: value for key, value in data.items() if key != "items"}
        items.extend(data.get("items", []))

    if merged is None:
        raise SystemExit("No manifest content loaded")

    merged["items"] = items
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(merged, indent=2) + "\n", encoding="utf-8")

    print(f"Merged {len(manifest_paths)} manifest(s) into {output_path}")
    return 0
