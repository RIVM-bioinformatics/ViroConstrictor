"""Basic tests for manifest creation, validation, and JSON I/O helpers."""

import json
from pathlib import Path

from container_manager.src.merge_manifests import merge_manifests
from container_manager.src.models import Manifest, ManifestItem
from container_manager.src.publish_to_ghcr import validate_publish_manifest
from container_manager.src.version import REPO_ROOT


def test_manifest_creation_to_dict() -> None:
    """Create a manifest model and ensure it serializes expected fields."""
    item = ManifestItem(
        name="Clean",
        hash="abc123",
        docker_tag="samplepkg_clean:abc123",
        docker_ref="ghcr.io/rivm-bioinformatics/samplepkg_clean:abc123",
        dockerfile_path="container_manager/Clean.dockerfile",
        env_file="workflow/envs/Clean.yaml",
        artifact_tar="container_manager/samplepkg_clean_abc123.tar",
        artifact_sif="container_manager/samplepkg_clean_abc123.sif",
        build_status="built",
    )
    manifest = Manifest(
        schema_version="1",
        generated_at_utc="2026-01-01T00:00:00+00:00",
        git_commit="deadbeef",
        registry="ghcr.io/rivm-bioinformatics",
        items=[item],
    )

    payload = manifest.to_dict()

    assert payload["registry"] == "ghcr.io/rivm-bioinformatics"
    assert payload["items"][0]["docker_tag"] == "samplepkg_clean:abc123"


def test_manifest_validation_publish_ok() -> None:
    """Validate a minimal publish manifest payload."""
    artifact_tar = REPO_ROOT / "container_manager" / "samplepkg_clean_abc123.tar"
    manifest = {
        "registry": "ghcr.io/rivm-bioinformatics",
        "items": [
            {
                "name": "Clean",
                "build_status": "built",
                "docker_tag": "samplepkg_clean:abc123",
                "docker_ref": "ghcr.io/rivm-bioinformatics/samplepkg_clean:abc123",
                "artifact_tar": str(artifact_tar),
            }
        ],
    }

    validated = validate_publish_manifest(manifest)

    assert validated["items"][0]["name"] == "Clean"


def test_merge_manifests_combines_items(tmp_path: Path) -> None:
    """Merge two manifest files and verify combined items in output."""
    manifest_a = tmp_path / "a.json"
    manifest_b = tmp_path / "b.json"
    merged_path = tmp_path / "merged" / "out.json"

    manifest_a.write_text(
        json.dumps(
            {
                "schema_version": "1",
                "generated_at_utc": "2026-01-01T00:00:00+00:00",
                "git_commit": "aaaa",
                "registry": "ghcr.io/rivm-bioinformatics",
                "items": [{"name": "Clean"}],
            }
        ),
        encoding="utf-8",
    )
    manifest_b.write_text(
        json.dumps(
            {
                "schema_version": "1",
                "generated_at_utc": "2026-01-01T00:00:00+00:00",
                "git_commit": "bbbb",
                "registry": "ghcr.io/rivm-bioinformatics",
                "items": [{"name": "Alignment"}],
            }
        ),
        encoding="utf-8",
    )

    rc = merge_manifests(str(tmp_path / "*.json"), merged_path)

    assert rc == 0
    assert merged_path.exists()
    payload = json.loads(merged_path.read_text(encoding="utf-8"))
    assert len(payload["items"]) == 2
    assert payload["items"][0]["name"] == "Clean"
    assert payload["items"][1]["name"] == "Alignment"
