"""Publish built container images to the configured remote registry from a manifest."""

import re
import subprocess
from pathlib import Path
from typing import Any, cast

from container_manager.src.io import ensure_safe_manifest_path, ensure_safe_path, read_json, write_json
from container_manager.src.logging import get_logger
from container_manager.src.version import REPO_ROOT

logger = get_logger(__name__)

CONTAINER_MANAGER_DIR = REPO_ROOT / "container_manager"

_SAFE_NAME_RE = re.compile(r"^[A-Za-z0-9._-]+$")
_SAFE_TAG_RE = re.compile(r"^[a-z0-9][a-z0-9._-]*:[A-Za-z0-9][A-Za-z0-9._-]{0,127}$")
_SAFE_REGISTRY_RE = re.compile(r"^[a-z0-9.-]+(?::\d+)?(?:/[a-z0-9._-]+)*$")
_SAFE_DOCKER_REF_RE = re.compile(r"^[a-z0-9.-]+(?::\d+)?/[a-z0-9._/-]+:[A-Za-z0-9][A-Za-z0-9._-]{0,127}$")


def _require_string(mapping: dict[str, Any], key: str, context: str) -> str:
    value = mapping.get(key)
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"Invalid manifest {context}: field '{key}' must be a non-empty string")
    return value


def _validate_registry(registry: str) -> str:
    reg = registry.strip().lower()
    if not _SAFE_REGISTRY_RE.fullmatch(reg):
        raise ValueError(f"Invalid manifest registry: {registry}")
    return reg


def _validate_publish_manifest_item(item: dict[str, Any], registry: str | None, idx: int) -> None:
    context = f"item[{idx}]"
    name = _require_string(item, "name", context)
    if not _SAFE_NAME_RE.fullmatch(name):
        raise ValueError(f"Invalid manifest {context}: unsupported item name '{name}'")

    build_status = _require_string(item, "build_status", context)
    if build_status in {"built", "skipped_existing"}:
        docker_tag = _require_string(item, "docker_tag", context)
        docker_ref = _require_string(item, "docker_ref", context)
        if not _SAFE_TAG_RE.fullmatch(docker_tag):
            raise ValueError(f"Invalid manifest {context}: malformed docker_tag '{docker_tag}'")
        if not _SAFE_DOCKER_REF_RE.fullmatch(docker_ref):
            raise ValueError(f"Invalid manifest {context}: malformed docker_ref '{docker_ref}'")
        if registry and docker_ref != f"{registry}/{docker_tag}":
            raise ValueError(
                f"Invalid manifest {context}: docker_ref '{docker_ref}' does not match registry '{registry}' and docker_tag '{docker_tag}'"
            )

        artifact_tar = item.get("artifact_tar")
        if artifact_tar is not None:
            if not isinstance(artifact_tar, str) or not artifact_tar.strip():
                raise ValueError(f"Invalid manifest {context}: artifact_tar must be a non-empty string when provided")
            ensure_safe_path(Path(artifact_tar), suffix=".tar", base_dir=CONTAINER_MANAGER_DIR, must_exist=False)


def validate_publish_manifest(manifest: dict[str, Any]) -> dict[str, Any]:
    """Validate manifest structure and publish-critical fields before docker operations."""
    raw_items = manifest.get("items")
    if not isinstance(raw_items, list):
        raise ValueError("Invalid manifest: top-level 'items' must be a list")
    items = cast(list[object], raw_items)

    registry = None
    raw_registry = manifest.get("registry")
    if raw_registry is not None:
        if not isinstance(raw_registry, str) or not raw_registry.strip():
            raise ValueError("Invalid manifest: top-level 'registry' must be a non-empty string when provided")
        registry = _validate_registry(raw_registry)

    for idx, item_obj in enumerate(items):
        if not isinstance(item_obj, dict):
            raise ValueError(f"Invalid manifest: item[{idx}] must be an object")
        item = cast(dict[str, Any], item_obj)
        _validate_publish_manifest_item(item, registry, idx)
    return manifest


def validate_convert_manifest(manifest: dict[str, Any]) -> dict[str, Any]:
    """Validate manifest structure and conversion-critical fields before apptainer operations."""
    raw_items = manifest.get("items")
    if not isinstance(raw_items, list):
        raise ValueError("Invalid manifest: top-level 'items' must be a list")
    items = cast(list[object], raw_items)

    for idx, item_obj in enumerate(items):
        if not isinstance(item_obj, dict):
            raise ValueError(f"Invalid manifest: item[{idx}] must be an object")
        item = cast(dict[str, Any], item_obj)
        context = f"item[{idx}]"
        build_status = item.get("build_status")
        if build_status in {"built", "skipped_existing"}:
            artifact_tar = _require_string(item, "artifact_tar", context)
            artifact_sif = _require_string(item, "artifact_sif", context)
            ensure_safe_path(Path(artifact_tar), suffix=".tar", base_dir=CONTAINER_MANAGER_DIR, must_exist=False)
            ensure_safe_path(Path(artifact_sif), suffix=".sif", base_dir=CONTAINER_MANAGER_DIR, must_exist=False)
    return manifest


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
    checked_manifest_path = ensure_safe_manifest_path(manifest_path, must_exist=True)
    logger.info("Starting publish from manifest: %s", checked_manifest_path)
    manifest = validate_publish_manifest(read_json(checked_manifest_path))
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

    write_json(checked_manifest_path, manifest)
    logger.info("Updated manifest after publish: %s", checked_manifest_path)
    return manifest
