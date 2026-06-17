"""Plan and build container artifacts, then record outcomes in a manifest."""

import os
import subprocess
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import cast

import requests  # type: ignore[import-untyped]

from container_manager.src.config import load_container_specs
from container_manager.src.hash import container_hash
from container_manager.src.io import ensure_safe_manifest_path, write_json
from container_manager.src.logging import get_logger
from container_manager.src.models import BuildPlanItem, Manifest, ManifestItem
from container_manager.src.version import HASH_SCHEMA_VERSION

logger = get_logger(__name__)


def _git_commit(repo_root: Path) -> str:
    """Return the current git commit hash for the repository, or 'unknown' if unavailable."""
    try:
        result = subprocess.run(
            ["git", "-C", str(repo_root), "rev-parse", "HEAD"],
            check=True,
            capture_output=True,
            text=True,
        )
        return result.stdout.strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        return "unknown"


def _existing_tags(container_name: str, token: str | None, ghcr_org: str | None) -> set[str]:
    """Fetch existing GHCR tags for a container package name using the GitHub API."""
    if not token or not ghcr_org:
        return set()
    endpoint = f"https://api.github.com/orgs/{ghcr_org}/packages/container/{container_name}/versions"
    headers = {
        "Accept": "application/vnd.github+json",
        "X-GitHub-Api-Version": "2022-11-28",
        "Authorization": f"Bearer {token}",
    }
    response = requests.get(endpoint, headers=headers, timeout=30)
    if response.status_code >= 400:
        return set()
    payload_obj: object = response.json()
    if not isinstance(payload_obj, list):
        return set()
    payload: list[object] = cast(list[object], payload_obj)
    tags = _process_payload(payload)
    return tags


def _process_payload(payload: list[object]) -> set[str]:
    tags: set[str] = set()
    for raw_version in payload:
        if not isinstance(raw_version, dict):
            continue
        version: dict[str, object] = cast(dict[str, object], raw_version)

        raw_metadata = version.get("metadata", {})
        if not isinstance(raw_metadata, dict):
            continue
        metadata: dict[str, object] = cast(dict[str, object], raw_metadata)

        raw_container = metadata.get("container", {})
        if not isinstance(raw_container, dict):
            continue
        container: dict[str, object] = cast(dict[str, object], raw_container)

        raw_tags = container.get("tags", [])
        if not isinstance(raw_tags, list):
            continue
        tags_list: list[object] = cast(list[object], raw_tags)
        for tag in tags_list:
            if isinstance(tag, str):
                tags.add(tag)
    return tags


def _registry_from_config(config: dict[str, object]) -> str:
    """Return configured registry value."""
    registry = config.get("registry")
    if not isinstance(registry, str) or not registry.strip():
        raise ValueError("Config field 'registry' must be a non-empty string")
    return registry.strip()


def _container_package_prefix_from_config(config: dict[str, object]) -> str:
    """Return configured package prefix used for image names."""
    package_prefix = config.get("container_package_prefix")
    if not isinstance(package_prefix, str) or not package_prefix.strip():
        raise ValueError("Config field 'container_package_prefix' must be a non-empty string")
    return package_prefix.strip().lower()


def _oci_labels_from_config(config: dict[str, object]) -> dict[str, str]:
    """Return OCI labels configured in YAML."""
    raw_labels = config.get("oci_labels")
    if raw_labels is None:
        raise ValueError("Config field 'oci_labels' is required and must be a mapping of string keys/values")
    if not isinstance(raw_labels, dict):
        raise ValueError("Config field 'oci_labels' must be a mapping of string keys/values")
    labels: dict[str, str] = {}
    for key, value in raw_labels.items():
        if not isinstance(key, str) or not isinstance(value, str):
            raise ValueError("Config field 'oci_labels' must contain string keys and string values")
        labels[key] = value
    return labels


def _ghcr_org_from_registry(registry: str) -> str | None:
    """Extract GHCR organization segment from registry string when present."""
    parts = registry.split("/", 2)
    if len(parts) >= 2 and parts[0].lower() == "ghcr.io" and parts[1]:
        return parts[1]
    return None


def plan_builds(
    repo_root: Path,
    config_path: Path,
    template_path: Path,
    output_dir: Path,
    only: set[str] | None = None,
) -> list[BuildPlanItem]:
    """Create a deterministic build plan for selected container specs."""
    config, specs = load_container_specs(config_path)
    registry = _registry_from_config(config)
    package_prefix = _container_package_prefix_from_config(config)
    plans: list[BuildPlanItem] = []
    for spec in specs:
        if only and spec.name not in only:
            continue
        env_path = repo_root / spec.env_file
        dockerfile_path = output_dir / spec.dockerfile_path
        hash_value = container_hash(spec=spec, env_path=env_path, template_path=template_path)
        docker_tag = f"{package_prefix}_{spec.name.lower()}:{hash_value}"
        docker_ref = f"{registry}/{docker_tag}"
        artifact_base = docker_tag.replace(":", "_")
        plans.append(
            BuildPlanItem(
                name=spec.name,
                hash=hash_value,
                docker_tag=docker_tag,
                docker_ref=docker_ref,
                env_file=str(env_path),
                dockerfile_path=str(dockerfile_path),
                artifact_tar=str(repo_root / "container_manager" / f"{artifact_base}.tar"),
                artifact_sif=str(repo_root / "container_manager" / f"{artifact_base}.sif"),
            )
        )
    return plans


def _append_labels(dockerfile_path: Path, hash_value: str, labels: dict[str, str]) -> str:
    """Append standard OCI metadata labels and version hash to a Dockerfile text."""
    text = dockerfile_path.read_text(encoding="utf-8")
    text += "\n"
    for key, value in labels.items():
        escaped_value = value.replace('"', '\\"')
        text += f'LABEL {key}="{escaped_value}"\n'
    text += f'LABEL version="{hash_value}"\n'
    return text


def _ensure_dockerfiles_exist(plans: list[BuildPlanItem]) -> None:
    """Validate that Dockerfiles have already been generated for all planned builds."""
    missing = [plan.dockerfile_path for plan in plans if not Path(plan.dockerfile_path).exists()]
    if not missing:
        return
    missing_list = ", ".join(missing)
    raise FileNotFoundError("Missing generated Dockerfile(s): " f"{missing_list}. " "Run the 'generate' command first.")


def build(
    repo_root: Path,
    config_path: Path,
    template_path: Path,
    output_dir: Path,
    manifest_path: Path,
    only: set[str] | None = None,
    dry_run: bool = False,
    skip_existing: bool = True,
) -> Manifest:
    """Build container tar artifacts from generated Dockerfiles and write a build manifest."""
    checked_manifest_path = ensure_safe_manifest_path(manifest_path, must_exist=False)
    logger.info("Starting build workflow")
    config, _ = load_container_specs(config_path)
    registry = _registry_from_config(config)
    labels = _oci_labels_from_config(config)
    ghcr_org = _ghcr_org_from_registry(registry)

    plans = plan_builds(repo_root=repo_root, config_path=config_path, template_path=template_path, output_dir=output_dir, only=only)
    if not dry_run:
        _ensure_dockerfiles_exist(plans)
    logger.info("Planned %d container build(s)", len(plans))

    token = os.environ.get("TOKEN")
    items: list[ManifestItem] = []
    for plan in plans:
        item = ManifestItem(
            name=plan.name,
            hash=plan.hash,
            docker_tag=plan.docker_tag,
            docker_ref=plan.docker_ref,
            dockerfile_path=plan.dockerfile_path,
            env_file=plan.env_file,
            artifact_tar=plan.artifact_tar,
            artifact_sif=plan.artifact_sif,
            build_status="planned" if dry_run else "pending",
        )
        if dry_run:
            items.append(item)
            continue

        container_pkg_name = plan.docker_tag.split(":", 1)[0]
        existing_tags: set[str] = _existing_tags(container_pkg_name, token, ghcr_org) if skip_existing else set()
        if plan.hash in existing_tags:
            item.build_status = "skipped_existing"
            logger.info("Skipping existing image tag: %s", plan.docker_tag)
            items.append(item)
            continue

        dockerfile_with_labels = _append_labels(Path(plan.dockerfile_path), plan.hash, labels)
        with tempfile.NamedTemporaryFile(mode="w", delete=False, encoding="utf-8") as temp_dockerfile:
            temp_dockerfile.write(dockerfile_with_labels)
            temp_path = temp_dockerfile.name

        subprocess.run(
            [
                "docker",
                "build",
                "-t",
                plan.docker_tag,
                "-f",
                temp_path,
                str(repo_root),
                "--network",
                "host",
                "--no-cache",
            ],
            check=True,
        )
        subprocess.run(["docker", "save", "-o", plan.artifact_tar, plan.docker_tag], check=True)
        subprocess.run(["docker", "rmi", plan.docker_tag], check=False)
        item.build_status = "built"
        logger.info("Built image artifact: %s", plan.artifact_tar)
        items.append(item)

    manifest = Manifest(
        schema_version=HASH_SCHEMA_VERSION,
        generated_at_utc=datetime.now(timezone.utc).isoformat(),
        git_commit=_git_commit(repo_root),
        registry=registry,
        items=items,
    )
    write_json(checked_manifest_path, manifest.to_dict())
    logger.info("Wrote build manifest: %s", checked_manifest_path)
    return manifest
