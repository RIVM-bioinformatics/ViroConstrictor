"""Basic tests for sync_local_cache helper behavior."""

import subprocess
from pathlib import Path

from container_manager.src import sync_local_cache as sync_module


def test_get_hash_returns_matching_value(monkeypatch) -> None:
    """get_hash should return the hash for a matching container name."""

    monkeypatch.setattr(sync_module, "_fetch_container_hashes", lambda _root: {"Clean": "abc123", "Alignment": "def456"})

    assert sync_module.get_hash("clean", repo_root=Path("/repo")) == "abc123"
    assert sync_module.get_hash("unknown", repo_root=Path("/repo")) is None


def test_download_containers_dry_run_skips_pull(monkeypatch, tmp_path: Path) -> None:
    """Dry-run download should not call subprocess pull and should return success."""

    def fail_subprocess(*_args, **_kwargs):
        raise AssertionError("subprocess.run should not be called in dry-run mode")

    monkeypatch.setattr(sync_module, "_containers_to_download", lambda _cache_dir, _repo_root: ["samplepkg_clean_abc123"])
    monkeypatch.setattr(sync_module, "_upstream_registry", lambda _repo_root: "ghcr.io/rivm-bioinformatics")
    monkeypatch.setattr(subprocess, "run", fail_subprocess)

    rc = sync_module.download_containers(tmp_path, dryrun=True)

    assert rc == 0
