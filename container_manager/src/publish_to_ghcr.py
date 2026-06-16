"""Publish built container images to the configured remote registry from a manifest."""

import subprocess
from pathlib import Path
from typing import Any

from container_manager.src.io import read_json, write_json
from container_manager.src.logging import get_logger

logger = get_logger(__name__)


def _ensure_image_loaded(docker_tag: str, artifact_tar: str | None) -> None:
    """Ensure a Docker image tag is present locally, loading it from tar if needed."""
    present = subprocess.run(["docker", "image", "inspect", docker_tag], check=False, capture_output=True)
    if present.returncode == 0:
        return
    if not artifact_tar:
        raise FileNotFoundError(f"Image {docker_tag} is not in local Docker and no tar artifact was provided")

    artifact_tar_path = Path(artifact_tar)
    if not artifact_tar_path.exists() and ":" in artifact_tar:
        fallback = Path(artifact_tar.replace(":", "_"))
        if fallback.exists():
            artifact_tar_path = fallback
    if not artifact_tar_path.exists():
        raise FileNotFoundError(f"Missing tar artifact for {docker_tag}: {artifact_tar}")

    subprocess.run(["docker", "image", "load", "-i", str(artifact_tar_path)], check=True)


def publish_from_manifest(manifest_path: Path, dry_run: bool = False, only: set[str] | None = None) -> dict[str, Any]:
    """Tag and push eligible built images listed in a manifest."""
    logger.info("Starting publish from manifest: %s", manifest_path)
    manifest = read_json(manifest_path)
    for item in manifest.get("items", []):
        name = item.get("name", "")
        if only and name not in only:
            continue

        build_status = item.get("build_status")
        if build_status not in {"built", "skipped_existing"}:
            if dry_run:
                item["publish_status"] = "skipped_unbuilt"
                logger.info("Skipping unbuilt image in dry-run: %s (build_status=%s)", name, build_status)
                continue
            raise ValueError(f"Cannot publish unbuilt item: {name}")

        if build_status == "skipped_existing":
            item["publish_status"] = "already_published"
            logger.info("Skipping already-published image: %s", item.get("docker_ref"))
            continue

        docker_tag = item["docker_tag"]
        docker_ref = item["docker_ref"]
        artifact_tar = item.get("artifact_tar")
        if dry_run:
            item["publish_status"] = "planned"
            continue

        _ensure_image_loaded(docker_tag, artifact_tar)
        subprocess.run(["docker", "tag", docker_tag, docker_ref], check=True)
        subprocess.run(["docker", "push", docker_ref], check=True)
        item["publish_status"] = "published"
        logger.info("Published image: %s", docker_ref)

    write_json(manifest_path, manifest)
    logger.info("Updated manifest after publish: %s", manifest_path)
    return manifest
