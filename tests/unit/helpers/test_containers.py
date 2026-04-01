"""Unit tests for :mod:`ViroConstrictor.workflow.helpers.containers`.

These tests focus on intended behavior with hermetic filesystem and subprocess
mocking. External tools (apptainer/singularity) are always mocked.
"""

from __future__ import annotations

import hashlib
from pathlib import Path

import pytest

import ViroConstrictor.workflow.helpers.containers as containers


def test_fetch_recipes_returns_absolute_yaml_paths(tmp_path: Path) -> None:
    """Return only YAML recipe files from one folder."""
    (tmp_path / "one.yaml").write_text("a: 1", encoding="utf-8")
    (tmp_path / "two.yaml").write_text("b: 2", encoding="utf-8")
    (tmp_path / "ignore.txt").write_text("x", encoding="utf-8")

    recipes = containers.fetch_recipes(str(tmp_path))

    assert len(recipes) == 2
    assert all(path.endswith(".yaml") for path in recipes)
    assert all(Path(path).is_absolute() for path in recipes)


def test_fetch_scripts_recurses_and_returns_python_files_only(tmp_path: Path) -> None:
    """Collect Python scripts recursively and ignore non-Python files."""
    nested = tmp_path / "nested"
    nested.mkdir()
    (tmp_path / "a.py").write_text("print('a')", encoding="utf-8")
    (nested / "b.py").write_text("print('b')", encoding="utf-8")
    (nested / "notes.md").write_text("# docs", encoding="utf-8")

    scripts = containers.fetch_scripts(str(tmp_path))

    assert len(scripts) == 2
    assert all(path.endswith(".py") for path in scripts)
    assert all(Path(path).is_absolute() for path in scripts)


def test_fetch_files_returns_all_entries_in_folder(tmp_path: Path) -> None:
    """Return all direct child entries in the given folder."""
    (tmp_path / "alpha.txt").write_text("a", encoding="utf-8")
    (tmp_path / "beta.yaml").write_text("b", encoding="utf-8")

    fetched = containers.fetch_files(str(tmp_path))

    assert {Path(path).name for path in fetched} == {"alpha.txt", "beta.yaml"}


def test_calculate_hashes_is_content_based_and_deterministic(tmp_path: Path) -> None:
    """Hash values should match SHA-256 first 6 characters of content."""
    file_a = tmp_path / "a.txt"
    file_b = tmp_path / "b.txt"
    file_a.write_text("same-content", encoding="utf-8")
    file_b.write_text("different-content", encoding="utf-8")

    result = containers.calculate_hashes([str(file_a), str(file_b)])

    assert result[str(file_a)] == hashlib.sha256(b"same-content").hexdigest()[:6]
    assert result[str(file_b)] == hashlib.sha256(b"different-content").hexdigest()[:6]
    assert result[str(file_a)] != result[str(file_b)]


def test_fetch_hashes_returns_sha_prefixes_for_all_recipe_files(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """fetch_hashes should compute 6-char SHA-256 prefixes per recipe file."""
    recipe_a = tmp_path / "A.yaml"
    recipe_b = tmp_path / "B.yaml"
    recipe_a.write_text("dependencies:\n  - python=3.11\n", encoding="utf-8")
    recipe_b.write_text("dependencies:\n  - pandas=2.3\n", encoding="utf-8")

    monkeypatch.setattr(containers, "fetch_recipes", lambda _folder: [str(recipe_a), str(recipe_b)])

    hashes = containers.fetch_hashes()

    assert hashes[str(recipe_a)] == hashlib.sha256(recipe_a.read_bytes()).hexdigest()[:6]
    assert hashes[str(recipe_b)] == hashlib.sha256(recipe_b.read_bytes()).hexdigest()[:6]


def test_log_containers_logs_dependency_versions_per_recipe(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Dependency details from YAML recipes should be emitted to debug logs."""
    recipe = tmp_path / "Clean.yaml"
    recipe.write_text("dependencies:\n  - python=3.11\n  - pip:\n      - pytest==8.0\n", encoding="utf-8")

    debug_messages: list[str] = []
    monkeypatch.setattr(containers.log, "debug", lambda message: debug_messages.append(message))

    containers.log_containers({str(recipe): "abc123"})

    assert len(debug_messages) == 1
    assert "recipe file" in debug_messages[0]
    assert "python=3.11" in debug_messages[0]
    assert "pip::pytest==8.0" in debug_messages[0]


def test_parse_dependencies_normalizes_supported_input_shapes() -> None:
    """Normalize plain strings, dict-based managers, and fallback objects."""
    dependencies = [
        "python=3.11",
        {"pip": ["pytest==8.0", "ruff==0.6"]},
        42,
    ]

    parsed = containers.parse_dependencies(dependencies)

    assert "python=3.11" in parsed
    assert "pip::pytest==8.0" in parsed
    assert "pip::ruff==0.6" in parsed
    assert "42" in parsed


def test_get_hash_returns_matching_hash_or_none(monkeypatch: pytest.MonkeyPatch) -> None:
    """Lookup should return a matching hash for container substring, else None."""
    monkeypatch.setattr(
        containers,
        "fetch_hashes",
        lambda: {
            "/tmp/Clean.yaml": "abc123",
            "/tmp/Consensus.yaml": "fff111",
        },
    )

    assert containers.get_hash("Clean") == "abc123"
    assert containers.get_hash("NotPresent") is None


def test_containerization_installed_prefers_apptainer_and_falls_back_singularity(monkeypatch: pytest.MonkeyPatch) -> None:
    """Installed check should pass when either executable is available."""

    def fake_system_apptainer_present(cmd: str) -> int:
        """Mock system call: return 0 for apptainer, 1 otherwise.

        Parameters
        ----------
        cmd : str
            Command string to check.

        Returns
        -------
        int
            Exit code: 0 if 'which apptainer', else 1.
        """
        return 0 if cmd == "which apptainer" else 1

    monkeypatch.setattr(containers.os, "system", fake_system_apptainer_present)
    assert containers.containerization_installed() is True

    calls: list[str] = []

    def fake_system_singularity_only(cmd: str) -> int:
        """Mock system call: apptainer fails, singularity succeeds.

        Parameters
        ----------
        cmd : str
            Command string to check.

        Returns
        -------
        int
            Exit code: 0 for 'which singularity', 1 otherwise.
        """
        calls.append(cmd)
        if cmd == "which apptainer":
            return 1
        if cmd == "which singularity":
            return 0
        return 1

    monkeypatch.setattr(containers.os, "system", fake_system_singularity_only)
    assert containers.containerization_installed() is True
    assert calls == ["which apptainer", "which singularity"]


def test_containerization_installed_returns_false_when_both_missing(monkeypatch: pytest.MonkeyPatch) -> None:
    """Installed check should return False when neither executable exists."""
    monkeypatch.setattr(containers.os, "system", lambda _cmd: 1)
    assert containers.containerization_installed() is False


def test_containers_to_download_returns_only_missing_and_creates_directory(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Required container set should be compared against local .sif cache names."""
    apptainer_dir = tmp_path / "cache"

    monkeypatch.setattr(
        containers,
        "fetch_hashes",
        lambda: {
            "/x/Alignment.yaml": "111aaa",
            "/x/Clean.yaml": "222bbb",
        },
    )
    monkeypatch.setattr(containers, "log_containers", lambda _hashes: None)

    apptainer_dir.mkdir()
    (apptainer_dir / "viroconstrictor_clean_222bbb.sif").write_text("", encoding="utf-8")

    to_download = containers.containers_to_download(apptainer_dir)

    assert to_download == ["viroconstrictor_alignment_111aaa"]


def test_containers_to_download_creates_missing_cache_dir(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Cache directory should be created when absent before listing entries."""
    apptainer_dir = tmp_path / "new-cache"

    monkeypatch.setattr(containers, "fetch_hashes", lambda: {"/x/Clean.yaml": "abcdef"})
    monkeypatch.setattr(containers, "log_containers", lambda _hashes: None)

    result = containers.containers_to_download(apptainer_dir)

    assert apptainer_dir.exists()
    assert result == ["viroconstrictor_clean_abcdef"]


def test_containerization_executable_selects_apptainer_when_available(monkeypatch: pytest.MonkeyPatch) -> None:
    """Apptainer should be selected when discovery command succeeds."""
    monkeypatch.setattr(containers.subprocess, "call", lambda *args, **kwargs: 0)
    assert containers.containerization_executable() == "apptainer"


def test_containerization_executable_selects_singularity_when_apptainer_missing(monkeypatch: pytest.MonkeyPatch) -> None:
    """Singularity should be selected when apptainer discovery fails."""
    monkeypatch.setattr(containers.subprocess, "call", lambda *args, **kwargs: 1)
    assert containers.containerization_executable() == "singularity"


@pytest.mark.xfail(reason="Intended behavior: raise a clear error if neither runtime is installed; current implementation silently falls back to singularity.")
def test_containerization_executable_raises_when_no_runtime_available(monkeypatch: pytest.MonkeyPatch) -> None:
    """Expect explicit failure when no container runtime is present (intended)."""
    monkeypatch.setattr(containers.subprocess, "call", lambda *args, **kwargs: 1)
    with pytest.raises(RuntimeError, match="No supported container runtime"):
        containers.containerization_executable()


def test_download_containers_dryrun_logs_targets_without_subprocess(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """Dry run should report translated image tags and not execute pulls."""
    monkeypatch.setattr(containers, "containers_to_download", lambda _path: ["viroconstrictor_clean_abc123"])

    info_messages: list[str] = []
    monkeypatch.setattr(containers.log, "info", lambda message: info_messages.append(message))

    called_subprocess = {"called": False}

    def fake_call(*_args, **_kwargs) -> int:
        """Mock subprocess call that records invocation.

        Returns
        -------
        int
            Always returns 0 to indicate success.
        """
        called_subprocess["called"] = True
        return 0

    monkeypatch.setattr(containers.subprocess, "call", fake_call)

    exit_code = containers.download_containers(tmp_path, dryrun=True)

    assert exit_code == 0
    assert any("viroconstrictor_clean:abc123" in message for message in info_messages)
    assert called_subprocess["called"] is False


def test_download_containers_success_pulls_each_required_image(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """Download should call pull command for each translated image tag."""
    monkeypatch.setattr(
        containers,
        "containers_to_download",
        lambda _path: ["viroconstrictor_clean_abc123", "viroconstrictor_alignment_def456"],
    )
    monkeypatch.setattr(containers, "containerization_executable", lambda: "apptainer")

    commands: list[str] = []

    def fake_call(command: str, shell: bool, stderr, stdout) -> int:
        """Mock subprocess call that records command invocation.

        Parameters
        ----------
        command : str
            Command string being executed.
        shell : bool
            Whether command is executed in shell mode.
        stderr
            Standard error stream (mocked).
        stdout
            Standard output stream (mocked).

        Returns
        -------
        int
            Always returns 0 to indicate success.
        """
        assert shell is True
        commands.append(command)
        return 0

    monkeypatch.setattr(containers.subprocess, "call", fake_call)

    exit_code = containers.download_containers(tmp_path, dryrun=False, verbose=False)

    assert exit_code == 0
    assert len(commands) == 2
    assert "docker://ghcr.io/rivm-bioinformatics/viroconstrictor_clean:abc123" in commands[0]
    assert "docker://ghcr.io/rivm-bioinformatics/viroconstrictor_alignment:def456" in commands[1]


def test_download_containers_returns_one_and_logs_error_on_pull_failure(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """A failed pull should stop iteration and return non-zero status."""
    monkeypatch.setattr(containers, "containers_to_download", lambda _path: ["viroconstrictor_clean_abc123"])
    monkeypatch.setattr(containers, "containerization_executable", lambda: "apptainer")
    monkeypatch.setattr(containers.subprocess, "call", lambda *args, **kwargs: 1)

    errors: list[str] = []
    monkeypatch.setattr(containers.log, "error", lambda message: errors.append(message))

    exit_code = containers.download_containers(tmp_path, dryrun=False)

    assert exit_code == 1
    assert any("Failed to download container" in msg for msg in errors)


def test_construct_container_bind_args_deduplicates_existing_parent_paths(tmp_path: Path) -> None:
    """Bind args should include unique directories for existing sample files only."""
    sample_dir = tmp_path / "sample"
    sample_dir.mkdir()
    existing_fastq = sample_dir / "reads.fastq.gz"
    existing_fastq.write_text("@r\nACGT\n+\n!!!!\n", encoding="utf-8")

    samples = {
        "s1": {
            "READS": str(existing_fastq),
            "MISSING": str(sample_dir / "missing.fastq.gz"),
        },
        "s2": {
            "READS2": str(existing_fastq),
        },
    }

    bind_args = containers.construct_container_bind_args(samples)

    assert f"--bind {sample_dir}" in bind_args
    assert bind_args.count(f"--bind {sample_dir}") == 1
