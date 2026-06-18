"""Define typed data models for container planning and manifest records."""

from dataclasses import asdict, dataclass, field
from typing import Any


@dataclass
class ContainerSpec:
    """Represent one container specification loaded from configuration."""

    name: str
    env_file: str
    dockerfile_path: str
    description: str
    extra_copies: list[dict[str, str]] = field(default_factory=list)


@dataclass
class BuildPlanItem:
    """Represent one computed build target with artifact paths and tags."""

    name: str
    hash: str
    docker_tag: str
    docker_ref: str
    env_file: str
    dockerfile_path: str
    artifact_tar: str
    artifact_sif: str


@dataclass
class ManifestItem:
    """Represent one manifest entry tracking build and publish status."""

    name: str
    hash: str
    docker_tag: str
    docker_ref: str
    dockerfile_path: str
    env_file: str
    artifact_tar: str
    artifact_sif: str | None = None
    build_status: str = "pending"
    publish_status: str = "pending"
    image_digest: str | None = None


@dataclass
class Manifest:
    """Represent the full manifest payload for a platform command run."""

    schema_version: str
    generated_at_utc: str
    git_commit: str
    registry: str
    items: list[ManifestItem]

    def to_dict(self) -> dict[str, Any]:
        """Serialize the manifest dataclass tree into a plain dictionary."""
        return asdict(self)
