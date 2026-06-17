"""Basic tests for build container planning and dry-run manifest creation."""

from pathlib import Path
from typing import Any

from container_manager.src import build_containers
from container_manager.src.models import BuildPlanItem, ContainerSpec


def test_plan_builds_creates_expected_item(monkeypatch) -> None:
    """Plan builds should produce deterministic tags and artifact paths."""

    def fake_load_container_specs(_config_path: Path) -> tuple[dict[str, Any], list[ContainerSpec]]:
        config = {
            "registry": "ghcr.io/rivm-bioinformatics",
            "container_package_prefix": "samplepkg",
        }
        specs = [
            ContainerSpec(
                name="Clean",
                env_file="workflow/envs/Clean.yaml",
                dockerfile_path="container_manager/Clean.dockerfile",
                description="Clean env",
            )
        ]
        return config, specs

    monkeypatch.setattr(build_containers, "load_container_specs", fake_load_container_specs)
    monkeypatch.setattr(build_containers, "container_hash", lambda **_kwargs: "abc123")

    plans = build_containers.plan_builds(
        repo_root=Path("/repo"),
        config_path=Path("/repo/container_manager/config_dockerfiles.yaml"),
        template_path=Path("/repo/container_manager/Dockerfile.j2"),
        dockerfiles_dir=Path("/repo/dockerfiles"),
        output_dir=Path("/repo"),
    )

    assert len(plans) == 1
    plan = plans[0]
    assert plan.docker_tag == "samplepkg_clean:abc123"
    assert plan.docker_ref == "ghcr.io/rivm-bioinformatics/samplepkg_clean:abc123"
    assert plan.artifact_tar == "/repo/samplepkg_clean_abc123.tar"


def test_build_dry_run_creates_planned_manifest(monkeypatch, tmp_path: Path) -> None:
    """Dry-run build should create a manifest payload with planned statuses."""
    captured: dict[str, object] = {}
    manifest_path = tmp_path / "manifest.json"

    def fake_load_container_specs(_config_path: Path) -> tuple[dict[str, object], list[ContainerSpec]]:
        return ({"registry": "ghcr.io/rivm-bioinformatics", "oci_labels": {}}, [])

    def fake_write_json(path: Path, payload: dict[str, object]) -> None:
        captured["path"] = path
        captured["payload"] = payload

    monkeypatch.setattr(build_containers, "load_container_specs", fake_load_container_specs)
    monkeypatch.setattr(
        build_containers,
        "plan_builds",
        lambda **_kwargs: [
            BuildPlanItem(
                name="Clean",
                hash="abc123",
                docker_tag="samplepkg_clean:abc123",
                docker_ref="ghcr.io/rivm-bioinformatics/samplepkg_clean:abc123",
                env_file="workflow/envs/Clean.yaml",
                dockerfile_path="container_manager/Clean.dockerfile",
                artifact_tar="/repo/container_manager/samplepkg_clean_abc123.tar",
                artifact_sif="/repo/container_manager/samplepkg_clean_abc123.sif",
            )
        ],
    )
    monkeypatch.setattr(build_containers, "ensure_safe_manifest_path", lambda _path, must_exist=False: manifest_path)
    monkeypatch.setattr(build_containers, "write_json", fake_write_json)
    monkeypatch.setattr(build_containers, "_git_commit", lambda _repo_root: "deadbeef")

    manifest = build_containers.build(
        repo_root=Path("/repo"),
        config_path=Path("/repo/container_manager/config_dockerfiles.yaml"),
        template_path=Path("/repo/container_manager/Dockerfile.j2"),
        dockerfiles_dir=Path("/repo/dockerfiles"),
        output_dir=Path("/repo"),
        manifest_path=Path("ignored/by/mock.json"),
        dry_run=True,
    )

    assert manifest.items[0].build_status == "planned"
    assert captured["path"] == manifest_path
    payload = captured["payload"]
    assert isinstance(payload, dict)
    assert payload["registry"] == "ghcr.io/rivm-bioinformatics"


def test_build_non_dry_run_fails_early_when_docker_unavailable(monkeypatch) -> None:
    """Non-dry-run build should fail early when Docker preflight check fails."""
    manifest_path = Path("/repo/container_manager/manifests/out.json")

    def fake_load_container_specs(_config_path: Path) -> tuple[dict[str, object], list[ContainerSpec]]:
        return ({"registry": "ghcr.io/rivm-bioinformatics", "oci_labels": {}}, [])

    monkeypatch.setattr(build_containers, "load_container_specs", fake_load_container_specs)
    monkeypatch.setattr(build_containers, "ensure_safe_manifest_path", lambda _path, must_exist=False: manifest_path)
    monkeypatch.setattr(
        build_containers,
        "_ensure_docker_available",
        lambda: (_ for _ in ()).throw(EnvironmentError("docker missing")),
    )

    try:
        build_containers.build(
            repo_root=Path("/repo"),
            config_path=Path("/repo/container_manager/config_dockerfiles.yaml"),
            template_path=Path("/repo/container_manager/Dockerfile.j2"),
            dockerfiles_dir=Path("/repo/dockerfiles"),
            output_dir=Path("/repo"),
            manifest_path=manifest_path,
            dry_run=False,
        )
    except EnvironmentError as exc:
        assert "docker missing" in str(exc)
    else:
        raise AssertionError("Expected EnvironmentError when Docker preflight fails")
