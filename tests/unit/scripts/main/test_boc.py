"""
Unit tests for Boc (breadth of coverage) script.

Tests calculation of breadth of coverage at multiple threshold levels
(1, 5, 10, 50, 100×) for single reference sequences.
"""

import sys
from argparse import ArgumentParser
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(PROJECT_ROOT.joinpath("ViroConstrictor/workflow")))
from ViroConstrictor.workflow.main.scripts.boc import Boc  # isort:skip


def _write_fasta(path: Path, records: list[tuple[str, str]]) -> Path:
    """
    Write FASTA format reference file to path.

    Parameters
    ----------
    path : Path
        Destination file path.
    records : list[tuple[str, str]]
        List of (sequence_id, sequence) tuples.

    Returns
    -------
    Path
        Path to the written FASTA file.
    """
    contents = "".join(f">{name}\n{sequence}\n" for name, sequence in records)
    path.write_text(contents, encoding="utf-8")
    return path


def _write_coverage(path: Path, depths: list[int]) -> Path:
    """
    Write per-position coverage file to path.

    Parameters
    ----------
    path : Path
        Destination file path.
    depths : list[int]
        List of coverage depth values (1-indexed positions).

    Returns
    -------
    Path
        Path to the written coverage file.
    """
    lines = ["position\tcoverage\n"]
    lines.extend(f"{idx}\t{depth}\n" for idx, depth in enumerate(depths, start=1))
    path.write_text("".join(lines), encoding="utf-8")
    return path


def _parse_result_line(output_path: Path) -> list[str]:
    """
    Parse single-line TSV result file.

    Parameters
    ----------
    output_path : Path
        Path to result file.

    Returns
    -------
    list[str]
        Tab-separated values from the result line.
    """
    line = output_path.read_text(encoding="utf-8").strip()
    return line.split("\t")


def test_add_arguments_parses_required_flags() -> None:
    """
    Verify that Boc CLI parser registers all required command-line flags.

    Tests that the argument parser accepts input, output, samplename, and coverage
    flags and correctly stores their values.
    """
    parser = ArgumentParser()
    Boc.add_arguments(parser)

    args = parser.parse_args(["--input", "ref.fasta", "--output", "boc.tsv", "--samplename", "sample-a", "--coverage", "cov.tsv"])

    assert args.input == "ref.fasta"
    assert args.output == "boc.tsv"
    assert args.samplename == "sample-a"
    assert args.coverage == "cov.tsv"


def test_run_writes_expected_breadth_percentages(tmp_path: Path) -> None:
    """
    Verify breadth of coverage calculation at standard thresholds.

    Tests that the script correctly calculates coverage breadth at 1×, 5×, 10×, 50×,
    and 100× depth thresholds for each position in the reference sequence.
    """
    fasta = _write_fasta(tmp_path / "ref.fasta", [("ref", "ACGTACGTAC")])
    coverage = _write_coverage(tmp_path / "coverage.tsv", [0, 1, 4, 5, 9, 10, 49, 50, 99, 100])
    output = tmp_path / "boc.tsv"

    Boc(input=fasta, output=output, samplename="sample-1", coverage=coverage).run()

    columns = _parse_result_line(output)
    assert columns[0] == "sample-1"
    assert columns[1:] == ["90.0", "70.0", "50.0", "30.0", "10.0"]


def test_run_outputs_zero_width_when_all_positions_are_below_thresholds(tmp_path: Path) -> None:
    """
    Verify that zero coverage positions produce 0% breadth at all thresholds.

    Tests that when all positions have zero depth, the breadth of coverage
    is reported as 0% for all depth thresholds.
    """
    fasta = _write_fasta(tmp_path / "ref.fasta", [("ref", "A" * 8)])
    coverage = _write_coverage(tmp_path / "coverage.tsv", [0, 0, 0, 0, 0, 0, 0, 0])
    output = tmp_path / "boc.tsv"

    Boc(input=fasta, output=output, samplename="sample-zero", coverage=coverage).run()

    columns = _parse_result_line(output)
    assert columns[0] == "sample-zero"
    assert columns[1:] == ["0.0", "0.0", "0.0", "0.0", "0.0"]


def test_run_raises_file_not_found_for_missing_coverage(tmp_path: Path) -> None:
    """
    Verify that missing coverage file raises FileNotFoundError.

    Tests that the script fails with FileNotFoundError when the coverage
    input file does not exist.
    """
    fasta = _write_fasta(tmp_path / "ref.fasta", [("ref", "ACGT")])
    output = tmp_path / "boc.tsv"

    with pytest.raises(FileNotFoundError):
        Boc(input=fasta, output=output, samplename="sample-missing", coverage=tmp_path / "missing.tsv").run()


@pytest.mark.xfail(reason="Implementation does not validate empty FASTA with a clear ValueError", strict=False)
def test_run_rejects_empty_fasta_with_clear_error(tmp_path: Path) -> None:
    """
    Test empty FASTA input rejection (XFAIL - intended behavior).

    Intended behavior: Empty FASTA input should fail fast with a clear validation
    error message, but the current implementation does not validate this clearly.
    """
    # Intended behavior: fail fast with a clear validation error for empty reference sequences.
    fasta = _write_fasta(tmp_path / "empty.fasta", [])
    coverage = _write_coverage(tmp_path / "coverage.tsv", [10, 10, 10])
    output = tmp_path / "boc.tsv"

    with pytest.raises(ValueError, match="empty|reference|FASTA"):
        Boc(input=fasta, output=output, samplename="sample-empty", coverage=coverage).run()


@pytest.mark.xfail(reason="Implementation raises IndexError instead of input validation error when coverage column is missing", strict=False)
def test_run_rejects_coverage_without_depth_column(tmp_path: Path) -> None:
    """
    Test malformed coverage input rejection (XFAIL - intended behavior).

    Intended behavior: Coverage file missing the depth column should raise ValueError,
    but the current implementation raises IndexError instead.
    """
    # Intended behavior: reject malformed coverage input with a clear ValueError.
    fasta = _write_fasta(tmp_path / "ref.fasta", [("ref", "ACGT")])
    coverage = tmp_path / "coverage.tsv"
    coverage.write_text("position\n1\n2\n3\n4\n", encoding="utf-8")
    output = tmp_path / "boc.tsv"

    with pytest.raises(ValueError, match="coverage|column|malformed"):
        Boc(input=fasta, output=output, samplename="sample-bad", coverage=coverage).run()
