"""Synchronize local container cache state with required remote container artifacts."""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

from container_manager.src.config import load_container_specs
from container_manager.src.hash import container_hash
from container_manager.src.logging import get_logger
from container_manager.src.version import REPO_ROOT

CONFIG_FILENAME = "config_dockerfiles.yaml"
logger = get_logger(__name__)


def _fetch_container_hashes(repo_root: Path) -> dict[str, str]:
    """Compute short hashes for all configured containers using platform hash rules."""
    config_path = repo_root / "container_manager" / CONFIG_FILENAME
    template_path = repo_root / "container_manager" / "Dockerfile.j2"
    _, specs = load_container_specs(config_path)

    hashes: dict[str, str] = {}
    for spec in specs:
        env_path = repo_root / spec.env_file
        hashes[spec.name] = container_hash(spec=spec, env_path=env_path, template_path=template_path)
    return hashes


def _container_package_prefix(repo_root: Path) -> str:
    """Return configured package prefix for container image names."""
    config_path = repo_root / "container_manager" / CONFIG_FILENAME
    config, _ = load_container_specs(config_path)
    prefix = config.get("container_package_prefix")
    if not isinstance(prefix, str) or not prefix.strip():
        raise ValueError("Config field 'container_package_prefix' must be a non-empty string")
    return prefix.strip().lower()


def _upstream_registry(repo_root: Path) -> str:
    """Return configured upstream registry for container pulls."""
    config_path = repo_root / "container_manager" / CONFIG_FILENAME
    config, _ = load_container_specs(config_path)
    registry = config.get("registry")
    if not isinstance(registry, str) or not registry.strip():
        raise ValueError("Config field 'registry' must be a non-empty string")
    return registry.strip()


def _required_container_ids(repo_root: Path) -> list[str]:
    """Return required local container ids formatted as name_hash strings."""
    recipe_hashes = _fetch_container_hashes(repo_root)
    package_prefix = _container_package_prefix(repo_root)
    required: list[str] = []
    for recipe_name, hash_value in recipe_hashes.items():
        required.append(f"{package_prefix}_{recipe_name}_{hash_value}".lower())
    return required


def get_hash(target_container: str, repo_root: Path | None = None) -> str | None:
    """Return the short hash for a workflow env file whose stem matches the given container name."""
    root = repo_root or Path(__file__).resolve().parents[2]
    recipe_hashes = _fetch_container_hashes(root)
    target = target_container.lower()
    return next((hash_value for recipe_name, hash_value in recipe_hashes.items() if target == recipe_name.lower()), None)


def _containers_to_download(cache_dir: Path, repo_root: Path) -> list[str]:
    """Determine which required containers are missing from the local cache."""
    cache_dir.mkdir(parents=True, exist_ok=True)
    required = _required_container_ids(repo_root)
    present = {item.stem for item in cache_dir.iterdir() if item.is_file()}
    return [container for container in required if container not in present]


def containerization_executable() -> str:
    """Return the preferred container executable name available on this host."""
    if shutil.which("apptainer"):
        return "apptainer"
    return "singularity"


def containerization_installed() -> bool:
    """Return True when either apptainer or singularity is available."""
    return shutil.which("apptainer") is not None or shutil.which("singularity") is not None


def download_containers(apptainer_path: Path | str, dryrun: bool = False) -> int:
    """Download missing workflow containers into the local cache directory."""
    apptainer_path = Path(apptainer_path)

    to_download = _containers_to_download(apptainer_path, REPO_ROOT)
    tags = [x.rsplit("_", 1)[0] + ":" + x.rsplit("_", 1)[1] for x in to_download]
    upstream_registry = _upstream_registry(REPO_ROOT)
    logger.info("Container sync resolved %d missing image(s)", len(tags))

    if dryrun:
        if tags:
            logger.info("Container(s) to download: %s", ", ".join(tags))
        else:
            logger.info("No containers need to be downloaded")
        return 0

    executable = containerization_executable()
    for container in tags:
        remote_ref = f"docker://{upstream_registry}/{container}"
        logger.info("Downloading container: %s", container)
        result = subprocess.run(
            [executable, "pull", "--dir", str(apptainer_path), remote_ref],
            check=False,
            capture_output=True,
            text=True,
        )
        if result.returncode != 0:
            stderr_msg = (result.stderr or "").strip().splitlines()
            stderr_preview = stderr_msg[-1] if stderr_msg else "no stderr output"
            logger.error(
                "Failed to download container from GHCR: %s. "
                "This usually means the required hash tag is not published yet. "
                "Expected local cache file stem: %s. "
                "Check that this tag exists in GHCR or build/publish matching containers first. "
                "Pull error: %s",
                remote_ref,
                container.replace(":", "_"),
                stderr_preview,
            )
            return 1
        logger.info("Successfully downloaded container: %s", container)

    return 0


def sync_cache(cache_dir: Path, dry_run: bool = False) -> int:
    """Download missing containers into a cache directory and return status code."""
    cache_dir.mkdir(parents=True, exist_ok=True)
    return download_containers(cache_dir, dryrun=dry_run)
