"""Basic tests for conversion flow from manifest to apptainer commands."""

import subprocess
from pathlib import Path

from container_manager.src import convert_docker_tar_to_apptainer_sif as convert_module
from container_manager.src.version import REPO_ROOT


def test_convert_from_manifest_dry_run_marks_planned(monkeypatch, tmp_path: Path) -> None:
    """Dry-run conversion should mark eligible items as planned and write the manifest."""
    checked_manifest_path = tmp_path / "manifest.json"
    captured: dict[str, object] = {}
    artifact_tar = REPO_ROOT / "container_manager" / "samplepkg_clean_abc123.tar"
    artifact_sif = REPO_ROOT / "container_manager" / "samplepkg_clean_abc123.sif"

    manifest = {
        "items": [
            {
                "name": "Clean",
                "build_status": "built",
                "artifact_tar": str(artifact_tar),
                "artifact_sif": str(artifact_sif),
            }
        ]
    }

    def fake_write_json(path: Path, payload: dict[str, object]) -> None:
        captured["path"] = path
        captured["payload"] = payload

    def fail_subprocess(*_args, **_kwargs):
        raise AssertionError("subprocess.run should not be called in dry-run mode")

    monkeypatch.setattr(convert_module, "ensure_safe_manifest_path", lambda _path, must_exist=True: checked_manifest_path)
    monkeypatch.setattr(convert_module, "read_json", lambda _path: manifest)
    monkeypatch.setattr(convert_module, "write_json", fake_write_json)
    monkeypatch.setattr(subprocess, "run", fail_subprocess)

    result = convert_module.convert_from_manifest(Path("ignored.json"), dry_run=True)

    assert result["items"][0]["convert_status"] == "planned"
    assert captured["path"] == checked_manifest_path


def test_convert_from_manifest_runs_apptainer_commands(monkeypatch, tmp_path: Path) -> None:
    """Non-dry-run conversion should execute expected apptainer commands."""
    checked_manifest_path = tmp_path / "manifest.json"
    commands: list[list[str]] = []

    manifest = {"items": [{"name": "Clean", "build_status": "built", "artifact_tar": "x.tar", "artifact_sif": "y.sif"}]}

    def fake_subprocess_run(cmd: list[str], check: bool) -> None:
        assert check is True
        commands.append(cmd)

    monkeypatch.setattr(convert_module, "ensure_safe_manifest_path", lambda _path, must_exist=True: checked_manifest_path)
    monkeypatch.setattr(convert_module, "read_json", lambda _path: manifest)
    monkeypatch.setattr(convert_module, "validate_convert_manifest", lambda payload: payload)
    monkeypatch.setattr(convert_module, "_check_manifest", lambda _item, _dry_run, _only: (Path("/tmp/input.tar"), "/tmp/output.sif"))
    monkeypatch.setattr(convert_module, "write_json", lambda _path, _payload: None)
    monkeypatch.setattr(subprocess, "run", fake_subprocess_run)

    result = convert_module.convert_from_manifest(Path("ignored.json"), dry_run=False)

    assert result["items"][0]["convert_status"] == "converted"
    assert commands == [
        ["apptainer", "cache", "clean", "-f"],
        ["apptainer", "build", "/tmp/output.sif", "docker-archive:///tmp/input.tar"],
    ]
