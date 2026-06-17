"""Compute deterministic hashes for files and container build specifications."""

import hashlib
import json
from pathlib import Path

from container_manager.src.constants import HASH_SCHEMA_VERSION
from container_manager.src.models import ContainerSpec


def _sha256_bytes(data: bytes) -> str:
    """Return the SHA-256 hexadecimal digest for raw bytes."""
    return hashlib.sha256(data).hexdigest()


def file_hash(path: Path) -> str:
    """Return the SHA-256 digest of a file's raw content."""
    return _sha256_bytes(path.read_bytes())


def container_hash(spec: ContainerSpec, env_path: Path, template_path: Path) -> str:
    """Return a short stable hash for a container based on all build-relevant inputs."""
    payload = {
        "schema": HASH_SCHEMA_VERSION,
        "env_hash": file_hash(env_path),
        "template_hash": file_hash(template_path),
        "spec": {
            "name": spec.name,
            "env_file": spec.env_file,
            "dockerfile_path": spec.dockerfile_path,
            "description": spec.description,
            "extra_copies": spec.extra_copies,
        },
    }
    digest = _sha256_bytes(json.dumps(payload, sort_keys=True).encode("utf-8"))
    return digest[:6]
