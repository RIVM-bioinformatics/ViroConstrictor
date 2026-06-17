"""Basic tests for publish manifest validation and publish flow."""

import subprocess
from pathlib import Path

from container_manager.src import publish_to_ghcr as publish_module
from container_manager.src.version import REPO_ROOT


def test_validate_publish_manifest_accepts_valid_payload() -> None:
    """Validation should accept a simple valid publish manifest."""
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

    validated = publish_module.validate_publish_manifest(manifest)

    assert validated["items"][0]["name"] == "Clean"


def test_publish_from_manifest_dry_run_marks_planned(monkeypatch, tmp_path: Path) -> None:
    """Dry-run publish should mark eligible items as planned and avoid docker calls."""
    checked_manifest_path = tmp_path / "manifest.json"
    captured: dict[str, object] = {}
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

    def fake_write_json(path: Path, payload: dict[str, object]) -> None:
        captured["path"] = path
        captured["payload"] = payload

    def fail_subprocess(*_args, **_kwargs):
        raise AssertionError("subprocess.run should not be called in dry-run mode")

    monkeypatch.setattr(publish_module, "ensure_safe_manifest_path", lambda _path, must_exist=True: checked_manifest_path)
    monkeypatch.setattr(publish_module, "read_json", lambda _path: manifest)
    monkeypatch.setattr(publish_module, "write_json", fake_write_json)
    monkeypatch.setattr(subprocess, "run", fail_subprocess)

    result = publish_module.publish_from_manifest(Path("ignored.json"), dry_run=True)

    assert result["items"][0]["publish_status"] == "planned"
    assert captured["path"] == checked_manifest_path
