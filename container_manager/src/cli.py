"""Expose the command-line interface for container platform workflows."""

import argparse
from pathlib import Path
from typing import Callable

from container_manager.src.build_containers import build
from container_manager.src.convert_docker_tar_to_apptainer_sif import convert_from_manifest
from container_manager.src.generate_dockerfiles import generate_dockerfiles
from container_manager.src.logging import get_logger
from container_manager.src.merge_manifests import merge_manifests
from container_manager.src.models import ManifestItem
from container_manager.src.publish_to_ghcr import publish_from_manifest
from container_manager.src.sync_local_cache import sync_cache
from container_manager.src.version import REPO_ROOT, VERSION

logger = get_logger(__name__)

DEFAULT_CONFIG = REPO_ROOT / "container_manager" / "config_dockerfiles.yaml"
DEFAULT_TEMPLATE = REPO_ROOT / "container_manager" / "Dockerfile.j2"
DEFAULT_OUTPUT_DIR = REPO_ROOT / "container_manager"
DEFAULT_MANIFEST = REPO_ROOT / "container_manager" / "manifests" / "container-build-manifest.json"
DEFAULT_MERGE_INPUT_GLOB = "./artifacts/container-build-manifest-*.json"


def _csv_to_set(raw: str | None) -> set[str] | None:
    """Convert a comma-separated string into a cleaned set of values."""
    if not raw:
        return None
    values = {part.strip() for part in raw.split(",") if part.strip()}
    return values or None


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="CLI for the platform package")
    parser.add_argument("--version", action="store_true", help="Show the version of the package")

    subparsers = parser.add_subparsers(dest="command")

    # generate
    generate_cmd = subparsers.add_parser(
        "generate",
        help="Generate Dockerfiles based on the configuration and template",
        description="Example: python -m container_manager generate --config path/to/config.yaml --template path/to/Dockerfile.j2 --output-dir path/to/output",
    )
    generate_cmd.add_argument("--config", type=Path, default=DEFAULT_CONFIG, help=f"Path to the configuration YAML file. Default: {DEFAULT_CONFIG}")
    generate_cmd.add_argument(
        "--template", type=Path, default=DEFAULT_TEMPLATE, help=f"Path to the Dockerfile Jinja2 template. Default: {DEFAULT_TEMPLATE}"
    )
    generate_cmd.add_argument(
        "--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR, help=f"Directory to output the generated Dockerfiles. Default: {DEFAULT_OUTPUT_DIR}"
    )
    generate_cmd.add_argument("--only", type=str, default=None, help="Comma-separated list of image names to process (e.g., 'image1,image2')")
    generate_cmd.add_argument("--dry-run", action="store_true", help="Print the generated Dockerfiles to stdout instead of writing to files")

    # build
    build_cmd = subparsers.add_parser(
        "build",
        help="Build container images based on the generated Dockerfiles and write a build manifest, which contains information about the built images.",
        description="Example: python -m container_manager build --config path/to/config.yaml --template path/to/Dockerfile.j2 --output-dir path/to/output --manifest path/to/manifest.json",
    )
    build_cmd.add_argument("--config", type=Path, default=DEFAULT_CONFIG, help=f"Path to the configuration YAML file. Default: {DEFAULT_CONFIG}")
    build_cmd.add_argument(
        "--template", type=Path, default=DEFAULT_TEMPLATE, help=f"Path to the Dockerfile Jinja2 template. Default: {DEFAULT_TEMPLATE}"
    )
    build_cmd.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=f"Directory to output the generated Dockerfiles and build artifacts. Default: {DEFAULT_OUTPUT_DIR}",
    )
    build_cmd.add_argument(
        "--manifest", type=Path, default=DEFAULT_MANIFEST, help=f"Path to write the build manifest JSON file. Default: {DEFAULT_MANIFEST}"
    )
    build_cmd.add_argument("--only", type=str, default=None, help="Comma-separated list of image names to build (e.g., 'image1,image2')")
    build_cmd.add_argument("--dry-run", action="store_true", help="Print the planned build actions to stdout instead of executing them")

    # convert
    convert_cmd = subparsers.add_parser("convert", help="Convert docker built artifacts into Singularity images based on the manifest")
    convert_cmd.add_argument(
        "--manifest", type=Path, default=DEFAULT_MANIFEST, help=f"Path to the build manifest JSON file. Default: {DEFAULT_MANIFEST}"
    )
    convert_cmd.add_argument("--only", type=str, default=None, help="Comma-separated list of image names to convert (e.g., 'image1,image2')")
    convert_cmd.add_argument("--dry-run", action="store_true", help="Print the planned conversion actions to stdout instead of executing them")

    # publish
    publish_cmd = subparsers.add_parser("publish", help="Publish images to the GHCR registry based on the manifest")
    publish_cmd.add_argument(
        "--manifest", type=Path, default=DEFAULT_MANIFEST, help=f"Path to the build manifest JSON file. Default: {DEFAULT_MANIFEST}"
    )
    publish_cmd.add_argument("--only", type=str, default=None, help="Comma-separated list of image names to publish (e.g., 'image1,image2')")
    publish_cmd.add_argument("--dry-run", action="store_true", help="Print the planned publish actions to stdout instead of executing them")

    # merge-manifests
    merge_cmd = subparsers.add_parser(
        "merge-manifests",
        help="Merge per-container manifest JSON files into one manifest",
    )
    merge_cmd.add_argument(
        "--input-glob",
        default=DEFAULT_MERGE_INPUT_GLOB,
        help=f"Glob pattern for input manifest files. Default: {DEFAULT_MERGE_INPUT_GLOB}",
    )
    merge_cmd.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_MANIFEST,
        help=f"Path for merged manifest output. Default: {DEFAULT_MANIFEST}",
    )

    # sync-cache
    sync_cmd = subparsers.add_parser(
        "sync-cache",
        help="Sync the local cache by computing required image hashes from the env file, Dockerfile template, and container specification; then pull any missing hash-tagged images from GHCR.",
    )
    sync_cmd.add_argument("--cache-dir", type=Path, required=True, help="Path to the cache directory.")
    sync_cmd.add_argument("--dry-run", action="store_true", help="Print the planned sync actions to stdout instead of executing them")

    # local
    local_cmd = subparsers.add_parser(
        "local",
        help="Run local operations using the specified cache and configuration. Runs generate > build > convert, sync",
    )
    local_cmd.add_argument(
        "--cache-dir", type=Path, required=True, help="Path to the cache directory. This will be used to store the built container images locally."
    )
    local_cmd.add_argument("--config", type=Path, default=DEFAULT_CONFIG, help=f"Path to the configuration file. Default: {DEFAULT_CONFIG}")
    local_cmd.add_argument("--template", type=Path, default=DEFAULT_TEMPLATE, help=f"Path to the template file. Default: {DEFAULT_TEMPLATE}")
    local_cmd.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR, help=f"Path to the output directory. Default: {DEFAULT_OUTPUT_DIR}")
    local_cmd.add_argument(
        "--manifest", type=Path, default=DEFAULT_MANIFEST, help=f"Path to the build manifest JSON file. Default: {DEFAULT_MANIFEST}"
    )
    local_cmd.add_argument("--only", type=str, default=None, help="Comma-separated list of image names to process (e.g., 'image1,image2')")
    local_cmd.add_argument("--dry-run", action="store_true", help="Print the planned local workflow actions instead of executing side effects")

    return parser.parse_args()


def _handle_generate(args: argparse.Namespace) -> int:
    """Handle the 'generate' command."""
    logger.info("Generating Dockerfiles")
    generated = generate_dockerfiles(
        config_path=args.config,
        template_path=args.template,
        output_dir=args.output_dir,
        dry_run=args.dry_run,
    )
    for path in generated:
        logger.info("%s", path)
    return 0


def _handle_build(args: argparse.Namespace) -> int:
    """Handle the 'build' command."""
    logger.info("Building container artifacts")
    manifest = build(
        repo_root=REPO_ROOT,
        config_path=args.config,
        template_path=args.template,
        output_dir=args.output_dir,
        manifest_path=args.manifest,
        only=_csv_to_set(args.only),
        dry_run=args.dry_run,
    )
    logger.info("Wrote manifest with %d item(s): %s", len(manifest.items), args.manifest)
    return 0


def _handle_convert(args: argparse.Namespace) -> int:
    """Handle the 'convert' command."""
    logger.info("Converting tar artifacts to SIF")
    convert_from_manifest(args.manifest, dry_run=args.dry_run, only=_csv_to_set(args.only))
    logger.info("Updated manifest: %s", args.manifest)
    return 0


def _handle_publish(args: argparse.Namespace) -> int:
    """Handle the 'publish' command."""
    logger.info("Publishing container images to GHCR")
    publish_from_manifest(args.manifest, dry_run=args.dry_run, only=_csv_to_set(args.only))
    logger.info("Updated manifest: %s", args.manifest)
    return 0


def _handle_merge_manifests(args: argparse.Namespace) -> int:
    """Handle the 'merge-manifests' command."""
    logger.info("Merging per-container manifests")
    return merge_manifests(input_glob=args.input_glob, output_path=args.output)


def _handle_sync_cache(args: argparse.Namespace) -> int:
    """Handle the 'sync-cache' command."""
    logger.info("Syncing local container cache")
    return sync_cache(cache_dir=args.cache_dir, dry_run=args.dry_run)


def _handle_local(args: argparse.Namespace) -> int:
    """Handle the 'local' command."""
    logger.info("Running local container workflow (build -> convert -> cache sync)")
    generate_dockerfiles(
        config_path=args.config,
        template_path=args.template,
        output_dir=args.output_dir,
        dry_run=args.dry_run,
    )
    manifest = build(
        repo_root=REPO_ROOT,
        config_path=args.config,
        template_path=args.template,
        output_dir=args.output_dir,
        manifest_path=args.manifest,
        only=_csv_to_set(args.only),
        dry_run=args.dry_run,
    )
    convert_from_manifest(args.manifest, dry_run=args.dry_run, only=_csv_to_set(args.only))
    logger.info("Wrote manifest with %d item(s): %s", len(manifest.items), args.manifest)

    if args.dry_run:
        logger.info("[dry-run] Would ensure cache directory exists: %s", args.cache_dir)
    else:
        args.cache_dir.mkdir(parents=True, exist_ok=True)

    for item in manifest.items:
        _move_sifs_to_cache(args, item)
    logger.info("Synced local cache: %s", args.cache_dir)

    return 0


def _move_sifs_to_cache(args: argparse.Namespace, item: ManifestItem) -> None:
    """
    Move the generated SIF file to the cache directory, and clean up the tar artifact if it exists.
    """
    if item.artifact_sif and Path(item.artifact_sif).exists():
        target = args.cache_dir / Path(item.artifact_sif).name
        if args.dry_run:
            logger.info("[dry-run] Would move SIF to cache: %s -> %s", item.artifact_sif, target)
        else:
            Path(item.artifact_sif).replace(target)
    if item.artifact_tar:
        artifact_tar = Path(item.artifact_tar)
        if artifact_tar.exists():
            if args.dry_run:
                logger.info("[dry-run] Would delete tar artifact: %s", artifact_tar)
            else:
                artifact_tar.unlink()


def generate_map() -> dict[str, Callable[[argparse.Namespace], int]]:
    """Map command names to their handler functions."""
    return {
        "generate": _handle_generate,
        "build": _handle_build,
        "convert": _handle_convert,
        "publish": _handle_publish,
        "merge-manifests": _handle_merge_manifests,
        "sync-cache": _handle_sync_cache,
        "local": _handle_local,
    }


def main() -> int:
    """CLI dispatcher."""
    args = parse_args()
    if args.version:
        logger.info("%s", VERSION)
        return 0

    command_map = generate_map()
    if args.command in command_map:
        handler = command_map[args.command]
        return handler(args)
    logger.error("No command provided. Use --help.")
    return 1
