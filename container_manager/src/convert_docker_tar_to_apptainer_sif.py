"""Convert built OCI tar artifacts into Apptainer SIF images from a manifest."""

import subprocess
from pathlib import Path
from typing import Any

from container_manager.src.io import read_json, write_json
from container_manager.src.logging import get_logger

logger = get_logger(__name__)


def convert_from_manifest(manifest_path: Path, dry_run: bool = False, only: set[str] | None = None) -> dict[str, Any]:
    """Convert eligible manifest items from docker-archive tar files to SIF files."""
    logger.info("Starting conversion from manifest: %s", manifest_path)
    manifest = read_json(manifest_path)
    items = manifest.get("items", [])
    for item in items:
        name = item.get("name", "")
        if only and name not in only:
            continue
        if item.get("build_status") not in {"built", "skipped_existing"}:
            continue
        artifact_tar = item.get("artifact_tar")
        artifact_sif = item.get("artifact_sif")
        if not artifact_tar or not artifact_sif:
            raise ValueError(f"Missing artifact paths for {name}")
        if dry_run:
            item["convert_status"] = "planned"
            continue
        artifact_tar_path = Path(artifact_tar)
        if not artifact_tar_path.exists() and ":" in artifact_tar:
            fallback = Path(artifact_tar.replace(":", "_"))
            if fallback.exists():
                artifact_tar_path = fallback
                item["artifact_tar"] = str(fallback)
        if not artifact_tar_path.exists():
            raise FileNotFoundError(f"Missing tar artifact: {artifact_tar}")
        subprocess.run(["apptainer", "cache", "clean", "-f"], check=True)
        subprocess.run(["apptainer", "build", artifact_sif, f"docker-archive://{artifact_tar_path}"], check=True)
        item["convert_status"] = "converted"
        logger.info("Converted artifact to SIF: %s", artifact_sif)

    write_json(manifest_path, manifest)
    logger.info("Updated manifest after conversion: %s", manifest_path)
    return manifest
