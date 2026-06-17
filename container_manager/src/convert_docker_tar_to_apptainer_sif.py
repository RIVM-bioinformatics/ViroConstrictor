"""Convert built OCI tar artifacts into Apptainer SIF images from a manifest."""

import shutil
import subprocess
from pathlib import Path
from typing import Any

from container_manager.src.io import ensure_safe_manifest_path, read_json, write_json
from container_manager.src.logging import get_logger
from container_manager.src.publish_to_ghcr import validate_convert_manifest

logger = get_logger(__name__)


def _ensure_apptainer_available() -> None:
    """Fail early when Apptainer is unavailable."""
    if shutil.which("apptainer") is None:
        raise EnvironmentError(
            "Apptainer executable was not found in PATH. Install Apptainer or run this step in an environment where "
            "the 'apptainer' command is available."
        )


def _check_manifest(manifest: dict[str, Any], dry_run: bool, only: set[str] | None) -> tuple[Path, str] | None:
    name = manifest.get("name", "unknown")
    if only and name not in only:
        return
    if manifest.get("build_status") not in {"built", "skipped_existing"}:
        return
    artifact_tar = manifest.get("artifact_tar")
    artifact_sif = manifest.get("artifact_sif")
    if not artifact_tar or not artifact_sif:
        raise ValueError(f"Missing artifact paths for {name}")
    if dry_run:
        manifest["convert_status"] = "planned"
        return
    artifact_tar_path = Path(artifact_tar)
    if not artifact_tar_path.exists() and ":" in artifact_tar:
        fallback = Path(artifact_tar.replace(":", "_"))
        if fallback.exists():
            artifact_tar_path = fallback
            manifest["artifact_tar"] = str(fallback)
    if not artifact_tar_path.exists():
        raise FileNotFoundError(f"Missing tar artifact: {artifact_tar}")
    return artifact_tar_path, artifact_sif


def convert_from_manifest(manifest_path: Path, dry_run: bool = False, only: set[str] | None = None) -> dict[str, Any]:
    """Convert eligible manifest items from docker-archive tar files to SIF files."""
    checked_manifest_path = ensure_safe_manifest_path(manifest_path, must_exist=True)
    logger.info("Starting conversion from manifest: %s", checked_manifest_path)
    manifest = validate_convert_manifest(read_json(checked_manifest_path))
    if not dry_run:
        _ensure_apptainer_available()
    items = manifest.get("items", [])
    for item in items:
        artifact_tar_path, artifact_sif = _check_manifest(item, dry_run, only) or (None, None)
        if artifact_tar_path is None or artifact_sif is None:
            continue
        subprocess.run(["apptainer", "cache", "clean", "-f"], check=True)
        subprocess.run(["apptainer", "build", artifact_sif, f"docker-archive://{artifact_tar_path}"], check=True)
        item["convert_status"] = "converted"
        logger.info("Converted artifact to SIF: %s", artifact_sif)

    write_json(checked_manifest_path, manifest)
    logger.info("Updated manifest after conversion: %s", checked_manifest_path)
    return manifest
