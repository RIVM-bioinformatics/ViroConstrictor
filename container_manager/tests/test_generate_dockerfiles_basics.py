"""Basic tests for Dockerfile generation from config and template."""

from pathlib import Path

from container_manager.src import generate_dockerfiles as generate_module
from container_manager.src.models import ContainerSpec


def test_generate_dockerfiles_dry_run_returns_paths(monkeypatch, tmp_path: Path) -> None:
    """Dry-run generation should return target paths without writing files."""

    def fake_load_container_specs(_config_path: Path) -> tuple[dict[str, object], list[ContainerSpec]]:
        config = {"base_image": "ubuntu:22.04"}
        specs = [
            ContainerSpec(
                name="Clean",
                env_file="workflow/envs/Clean.yaml",
                dockerfile_path="docker/Clean.dockerfile",
                description="Clean env",
            )
        ]
        return config, specs

    template_path = tmp_path / "Dockerfile.j2"
    template_path.write_text("FROM {{ base_image }}\n# {{ dockerfile.description }}\n", encoding="utf-8")
    output_dir = tmp_path / "out"

    monkeypatch.setattr(generate_module, "load_container_specs", fake_load_container_specs)

    rendered_paths = generate_module.generate_dockerfiles(
        config_path=tmp_path / "config.yaml",
        template_path=template_path,
        output_dir=output_dir,
        dry_run=True,
    )

    assert rendered_paths == [output_dir / "docker/Clean.dockerfile"]
    assert not (output_dir / "docker/Clean.dockerfile").exists()


def test_generate_dockerfiles_writes_rendered_file(monkeypatch, tmp_path: Path) -> None:
    """Generation should render and write Dockerfile content."""

    def fake_load_container_specs(_config_path: Path) -> tuple[dict[str, object], list[ContainerSpec]]:
        config = {"base_image": "ubuntu:22.04"}
        specs = [
            ContainerSpec(
                name="Clean",
                env_file="workflow/envs/Clean.yaml",
                dockerfile_path="docker/Clean.dockerfile",
                description="Clean env",
            )
        ]
        return config, specs

    template_path = tmp_path / "Dockerfile.j2"
    template_path.write_text(
        "FROM {{ base_image }}\nCOPY {{ dockerfile.env_file }} /env.yml\n# {{ dockerfile.description }}\n",
        encoding="utf-8",
    )
    output_dir = tmp_path / "out"

    monkeypatch.setattr(generate_module, "load_container_specs", fake_load_container_specs)

    rendered_paths = generate_module.generate_dockerfiles(
        config_path=tmp_path / "config.yaml",
        template_path=template_path,
        output_dir=output_dir,
        dry_run=False,
    )

    target = output_dir / "docker/Clean.dockerfile"
    assert rendered_paths == [target]
    assert target.exists()
    content = target.read_text(encoding="utf-8")
    assert "FROM ubuntu:22.04" in content
    assert "COPY workflow/envs/Clean.yaml /env.yml" in content
