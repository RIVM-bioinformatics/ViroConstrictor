"""Unit tests for workflow-only container helper utilities."""

from __future__ import annotations

from pathlib import Path

import ViroConstrictor.workflow.helpers.containers as containers


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


def test_construct_container_bind_args_ignores_non_string_values(tmp_path: Path) -> None:
    """Non-string sample values should not be converted into bind mount paths."""
    sample_dir = tmp_path / "sample"
    sample_dir.mkdir()
    existing_fastq = sample_dir / "reads.fastq.gz"
    existing_fastq.write_text("@r\nACGT\n+\n!!!!\n", encoding="utf-8")

    samples = {
        "s1": {
            "READS": str(existing_fastq),
            "COUNT": 123,
            "FLAGS": ["A", "B"],
        }
    }

    bind_args = containers.construct_container_bind_args(samples)

    assert f"--bind {sample_dir}" in bind_args
    assert "123" not in bind_args
    assert "['A', 'B']" not in bind_args
