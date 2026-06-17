"""Basic tests for Dockerfile config merge/validation and rendering."""

from container_manager.src.config import load_container_specs
from container_manager.src.constants import REPO_ROOT
from container_manager.src.generate_dockerfiles import generate_dockerfiles, join_config_yamls


def _normalized_dockerfile(text: str) -> str:
    """Normalize known non-functional template differences for fixture comparison."""
    lines = text.splitlines()
    filtered: list[str] = []
    for line in lines:
        stripped = line.strip()
        if stripped == "{}":
            continue
        if "micromamba install" in stripped or "micromamba clean" in stripped:
            continue
        filtered.append(line.rstrip())
    return "\n".join(filtered).strip() + "\n"


def test_yaml_merge_and_validation_with_real_configs() -> None:
    """Load real base/user YAML files and ensure merged config validates successfully."""
    config_path = REPO_ROOT / "container_manager" / "config_dockerfiles.yaml"

    merged = join_config_yamls(config_path)
    config, specs = load_container_specs(merged)

    assert config["base_image"] == "mambaorg/micromamba:latest"
    assert config["container_package_prefix"] == "viroconstrictor"
    assert len(specs) == 6
    assert any(spec.dockerfile_path.endswith("Clean.dockerfile") for spec in specs)


def test_load_container_specs_applies_dockerfile_entry_defaults() -> None:
    """Missing/empty dockerfile entry fields should be filled with defaults."""
    config = {
        "registry": "ghcr.io/example",
        "container_package_prefix": "samplepkg",
        "oci_labels": {"k": "v"},
        "dockerfiles": [
            {
                "env_file": "./workflow/envs/Alignment.yaml",
                "output": "",
                "description": "",
                "extra_copies": None,
            }
        ],
    }

    _, specs = load_container_specs(config)
    spec = specs[0]

    assert spec.dockerfile_path == "dockerfiles/Alignment.dockerfile"
    assert spec.env_file == "./workflow/envs/Alignment.yaml"
    assert spec.description == "Container image for Alignment."
    assert spec.extra_copies == []


def test_generate_dockerfiles_matches_real_test_data() -> None:
    """Generate Dockerfiles from real template/config and compare against test_data fixtures."""
    config_path = REPO_ROOT / "container_manager" / "config_dockerfiles.yaml"
    template_path = REPO_ROOT / "container_manager" / "docker_template_specific.j2"
    expected_dir = REPO_ROOT / "container_manager" / "tests" / "test_data"

    output_dir = REPO_ROOT / "container_manager" / "tests" / ".tmp_generated_dockerfiles"
    rendered_paths = generate_dockerfiles(config_path=config_path, template_path=template_path, output_dir=output_dir, dry_run=False)

    try:
        assert len(rendered_paths) == 6
        for generated_path in rendered_paths:
            fixture_path = expected_dir / generated_path.name
            assert fixture_path.exists(), f"Missing fixture: {fixture_path}"

            generated_text = generated_path.read_text(encoding="utf-8")
            fixture_text = fixture_path.read_text(encoding="utf-8")
            assert _normalized_dockerfile(generated_text) == _normalized_dockerfile(fixture_text)
    finally:
        for path in sorted(output_dir.rglob("*"), reverse=True):
            if path.is_file():
                path.unlink()
            elif path.is_dir():
                path.rmdir()
        if output_dir.exists():
            output_dir.rmdir()
