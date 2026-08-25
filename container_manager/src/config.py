"""Load and validate container platform YAML configuration inputs."""

from pathlib import Path
from typing import Any, cast

from container_manager.src.constants import (
    CONTAINER_MANAGER_BASE_DOCKERFILE_CONFIG_PATH,
    CONTAINER_MANAGER_DOCKERFILES_DIR,
    WORKFLOW_ENVS_DIR,
)
from container_manager.src.io import load_yaml
from container_manager.src.models import ContainerSpec


def join_config_yamls(config_path: Path) -> dict[str, Any]:
    """Load base config and override with user config, returning merged result."""
    base_config = load_yaml(CONTAINER_MANAGER_BASE_DOCKERFILE_CONFIG_PATH) if CONTAINER_MANAGER_BASE_DOCKERFILE_CONFIG_PATH.exists() else {}
    user_config = load_yaml(config_path)
    merged = {**base_config, **user_config}
    return merged


def load_container_specs(config_input: Path | dict[str, Any]) -> tuple[dict[str, Any], list[ContainerSpec]]:
    """Parse dockerfile entries from a path or mapping into typed container specifications."""
    if isinstance(config_input, Path):
        config = join_config_yamls(config_input)
    else:
        config = config_input

    registry = config.get("registry")
    if not isinstance(registry, str) or not registry.strip():
        raise ValueError("Config field 'registry' must be a non-empty string")

    package_prefix = config.get("container_package_prefix")
    if not isinstance(package_prefix, str) or not package_prefix.strip():
        raise ValueError("Config field 'container_package_prefix' must be a non-empty string")

    oci_labels_obj = config.get("oci_labels")
    if not isinstance(oci_labels_obj, dict):
        raise ValueError("Config field 'oci_labels' must be a mapping of string keys/values")
    oci_labels = cast(dict[object, object], oci_labels_obj)
    for key, value in oci_labels.items():
        if not isinstance(key, str) or not isinstance(value, str):
            raise ValueError("Config field 'oci_labels' must contain string keys and string values")

    entries_obj = config.get("dockerfiles")
    if not isinstance(entries_obj, list):
        raise ValueError("Config is missing 'dockerfiles' list")
    entries = cast(list[dict[str, Any]], entries_obj)

    specs = _generate_container_specs(entries)

    return config, specs


def _generate_container_specs(entries: list[dict[str, Any]]) -> list[ContainerSpec]:
    specs: list[ContainerSpec] = []
    for entry in entries:
        output = _default_output(entry.get("output"), entry.get("env_file"))
        env_file = _default_env_file(entry.get("env_file"), output)
        description = _default_description(entry.get("description"), output)
        extra_copies = _default_extra_copies(entry.get("extra_copies"))

        if not output.endswith(".dockerfile"):
            raise ValueError("dockerfiles.output must be a .dockerfile string")

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

    return specs


def _default_output(output: Any, env_file: Any) -> str:
    if isinstance(output, str) and output.strip():
        return output.strip()
    if isinstance(env_file, str) and env_file.strip():
        stem = Path(env_file.strip()).stem
        return str(CONTAINER_MANAGER_DOCKERFILES_DIR / f"{stem}.dockerfile")
    raise ValueError("dockerfiles.output must be provided, or dockerfiles.env_file must be set so output can be derived")


def _default_env_file(env_file: Any, output: str) -> str:
    if isinstance(env_file, str) and env_file.strip():
        return env_file.strip()
    stem = Path(output).stem
    return f"{WORKFLOW_ENVS_DIR}/{stem}.yaml"


def _default_description(description: Any, output: str) -> str:
    if isinstance(description, str) and description.strip():
        return description.strip()
    stem = Path(output).stem
    return f"Container image for {stem}."


def _default_extra_copies(extra_copies: Any) -> list[dict[str, str]]:
    if extra_copies is None:
        return []
    if isinstance(extra_copies, list):
        return cast(list[dict[str, str]], extra_copies)
    raise ValueError("dockerfiles.extra_copies must be a list")
