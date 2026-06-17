"""Render Dockerfiles from template and config entries for each container spec."""

from pathlib import Path
from typing import Any

from jinja2 import Environment, FileSystemLoader, StrictUndefined, select_autoescape

from container_manager.src.config import join_config_yamls, load_container_specs
from container_manager.src.io import ensure_dir
from container_manager.src.logging import get_logger

logger = get_logger(__name__)


def generate_dockerfiles(config_path: Path, template_path: Path, output_dir: Path, dry_run: bool = False) -> list[Path]:
    """Generate Dockerfiles for all configured container specs and return target paths."""
    logger.info("Generating Dockerfiles from %s", config_path)
    validated_raw_config = join_config_yamls(config_path)
    config, specs = load_container_specs(validated_raw_config)
    shared_fields = {key: value for key, value in config.items() if key != "dockerfiles"}

    env = Environment(
        loader=FileSystemLoader(str(template_path.parent)),
        # Dockerfiles are not HTML; only enable autoescaping for HTML/XML templates.
        # Effectively disabling it and ensuring that sonarqube doesnt complain about it.
        autoescape=select_autoescape(enabled_extensions=("html", "xml"), default_for_string=False, default=False),
        undefined=StrictUndefined,
        trim_blocks=True,
        lstrip_blocks=True,
    )
    template = env.get_template(template_path.name)

    rendered_paths: list[Path] = []
    ensure_dir(output_dir)
    for spec in specs:
        dockerfile_cfg: dict[str, Any] = {
            "output": spec.dockerfile_path,
            "env_file": spec.env_file,
            "description": spec.description,
            "extra_copies": spec.extra_copies,
        }
        rendered = template.render(**shared_fields, dockerfile=dockerfile_cfg).strip() + "\n"
        target = output_dir / spec.dockerfile_path
        rendered_paths.append(target)
        if dry_run:
            continue
        ensure_dir(target.parent)
        target.write_text(rendered, encoding="utf-8")

    logger.info("Prepared %d Dockerfile(s)", len(rendered_paths))

    return rendered_paths
