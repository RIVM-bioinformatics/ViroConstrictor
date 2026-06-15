import argparse
from pathlib import Path
from typing import Any

import yaml
from jinja2 import Environment, FileSystemLoader, StrictUndefined


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for Dockerfile generation."""
    parser = argparse.ArgumentParser(description="Render Dockerfiles from Jinja2 template and YAML config.")
    parser.add_argument(
        "--config",
        type=Path,
        default=Path("containers/dockerfiles.yaml"),
        help="Path to the YAML configuration file.",
    )
    parser.add_argument(
        "--template",
        type=Path,
        default=Path("containers/Dockerfile.j2"),
        help="Path to the Jinja2 Dockerfile template.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("containers"),
        help="Directory where rendered Dockerfiles are written.",
    )
    return parser.parse_args()


def load_config(config_path: Path) -> dict[str, Any]:
    """Load and validate YAML configuration.

    Parameters
    ----------
    config_path
        Path to the YAML configuration file.

    Returns
    -------
    dict[str, Any]
        Parsed configuration dictionary.
    """
    with config_path.open("r", encoding="utf-8") as handle:
        data: dict[str, Any] = yaml.safe_load(handle) or {}

    if "dockerfiles" not in data or not isinstance(data["dockerfiles"], list):
        raise ValueError("Configuration must include a 'dockerfiles' list.")

    return data


def render_dockerfiles(config: dict[str, Any], template_path: Path, output_dir: Path) -> list[Path]:
    """Render Dockerfiles from template and config.

    Parameters
    ----------
    config
        Parsed configuration dictionary.
    template_path
        Path to the Jinja2 template file.
    output_dir
        Output directory where Dockerfiles are written.

    Returns
    -------
    list[Path]
        Paths of rendered Dockerfiles.
    """
    environment = Environment(  # NOSONAR (S5247)
        loader=FileSystemLoader(str(template_path.parent)),
        undefined=StrictUndefined,
        trim_blocks=True,
        lstrip_blocks=True,
    )
    template = environment.get_template(template_path.name)

    shared_fields = {key: value for key, value in config.items() if key != "dockerfiles"}

    rendered_files: list[Path] = []
    output_dir.mkdir(parents=True, exist_ok=True)

    for dockerfile_cfg in config["dockerfiles"]:
        rendered = template.render(
            **shared_fields,
            dockerfile=dockerfile_cfg,
        ).strip()
        output_path = output_dir / dockerfile_cfg["output"]
        output_path.write_text(f"{rendered}\n", encoding="utf-8")
        rendered_files.append(output_path)

    return rendered_files


def main() -> None:
    """Generate Dockerfiles defined in the YAML configuration."""
    args = parse_args()
    config = load_config(args.config)
    rendered_files = render_dockerfiles(config=config, template_path=args.template, output_dir=args.output_dir)
    for path in rendered_files:
        print(f"Rendered {path}")


if __name__ == "__main__":
    main()
