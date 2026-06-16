"""Load and validate container platform YAML configuration inputs."""

from pathlib import Path
from typing import Any

from container_manager.src.io import load_yaml
from container_manager.src.models import ContainerSpec


def load_container_specs(config_path: Path) -> tuple[dict[str, Any], list[ContainerSpec]]:
    """Parse dockerfile entries from config into typed container specifications."""
    config = load_yaml(config_path)

    registry = config.get("registry")
    if not isinstance(registry, str) or not registry.strip():
        raise ValueError("Config field 'registry' must be a non-empty string")

    package_prefix = config.get("container_package_prefix")
    if not isinstance(package_prefix, str) or not package_prefix.strip():
        raise ValueError("Config field 'container_package_prefix' must be a non-empty string")

    oci_labels = config.get("oci_labels")
    if not isinstance(oci_labels, dict):
        raise ValueError("Config field 'oci_labels' must be a mapping of string keys/values")
    for key, value in oci_labels.items():
        if not isinstance(key, str) or not isinstance(value, str):
            raise ValueError("Config field 'oci_labels' must contain string keys and string values")

    entries: list[dict[str, str | dict[str, str]]] = config.get("dockerfiles")
    if not isinstance(entries, list):
        raise ValueError("Config is missing 'dockerfiles' list")

    specs: list[ContainerSpec] = []
    for entry in entries:
        output = entry.get("output")
        env_file = entry.get("env_file")
        description = entry.get("description", "")
        extra_copies = entry.get("extra_copies", [])
        if not isinstance(output, str) or not output.endswith(".dockerfile"):
            raise ValueError("dockerfiles.output must be a .dockerfile string")
        if not isinstance(env_file, str):
            raise ValueError("dockerfiles.env_file must be a string")
        if not isinstance(description, str):
            raise ValueError("dockerfiles.description must be a string")
        if not isinstance(extra_copies, list):
            raise ValueError("dockerfiles.extra_copies must be a list")

        name = Path(output).stem
        specs.append(
            ContainerSpec(
                name=name,
                env_file=env_file,
                dockerfile_path=output,
                description=description,
                extra_copies=extra_copies,
            )
        )

    return config, specs
