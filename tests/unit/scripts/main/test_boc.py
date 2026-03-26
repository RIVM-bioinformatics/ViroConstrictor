import sys
from argparse import ArgumentParser
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(PROJECT_ROOT.joinpath("ViroConstrictor/workflow")))
from ViroConstrictor.workflow.main.scripts.boc import Boc  # isort:skip


def _write_fasta(path: Path, records: list[tuple[str, str]]) -> Path:
    contents = "".join(f">{name}\n{sequence}\n" for name, sequence in records)
    path.write_text(contents, encoding="utf-8")
    return path


def _write_coverage(path: Path, depths: list[int]) -> Path:
    lines = ["position\tcoverage\n"]
    lines.extend(f"{idx}\t{depth}\n" for idx, depth in enumerate(depths, start=1))
    path.write_text("".join(lines), encoding="utf-8")
    return path


def _parse_result_line(output_path: Path) -> list[str]:
    line = output_path.read_text(encoding="utf-8").strip()
    return line.split("\t")


def test_add_arguments_parses_required_flags() -> None:
    parser = ArgumentParser()
    Boc.add_arguments(parser)

    args = parser.parse_args(["--input", "ref.fasta", "--output", "boc.tsv", "--samplename", "sample-a", "--coverage", "cov.tsv"])

    assert args.input == "ref.fasta"
    assert args.output == "boc.tsv"
    assert args.samplename == "sample-a"
    assert args.coverage == "cov.tsv"


def test_run_writes_expected_breadth_percentages(tmp_path: Path) -> None:
    fasta = _write_fasta(tmp_path / "ref.fasta", [("ref", "ACGTACGTAC")])
    coverage = _write_coverage(tmp_path / "coverage.tsv", [0, 1, 4, 5, 9, 10, 49, 50, 99, 100])
    output = tmp_path / "boc.tsv"

    Boc(input=fasta, output=output, samplename="sample-1", coverage=coverage).run()

    columns = _parse_result_line(output)
    assert columns[0] == "sample-1"
    assert columns[1:] == ["90.0", "70.0", "50.0", "30.0", "10.0"]


def test_run_outputs_zero_width_when_all_positions_are_below_thresholds(tmp_path: Path) -> None:
    fasta = _write_fasta(tmp_path / "ref.fasta", [("ref", "A" * 8)])
    coverage = _write_coverage(tmp_path / "coverage.tsv", [0, 0, 0, 0, 0, 0, 0, 0])
    output = tmp_path / "boc.tsv"

    Boc(input=fasta, output=output, samplename="sample-zero", coverage=coverage).run()

    columns = _parse_result_line(output)
    assert columns[0] == "sample-zero"
    assert columns[1:] == ["0.0", "0.0", "0.0", "0.0", "0.0"]


def test_run_raises_file_not_found_for_missing_coverage(tmp_path: Path) -> None:
    fasta = _write_fasta(tmp_path / "ref.fasta", [("ref", "ACGT")])
    output = tmp_path / "boc.tsv"

    with pytest.raises(FileNotFoundError):
        Boc(input=fasta, output=output, samplename="sample-missing", coverage=tmp_path / "missing.tsv").run()


@pytest.mark.xfail(reason="Implementation does not validate empty FASTA with a clear ValueError", strict=False)
def test_run_rejects_empty_fasta_with_clear_error(tmp_path: Path) -> None:
    # Intended behavior: fail fast with a clear validation error for empty reference sequences.
    fasta = _write_fasta(tmp_path / "empty.fasta", [])
    coverage = _write_coverage(tmp_path / "coverage.tsv", [10, 10, 10])
    output = tmp_path / "boc.tsv"

    with pytest.raises(ValueError, match="empty|reference|FASTA"):
        Boc(input=fasta, output=output, samplename="sample-empty", coverage=coverage).run()


@pytest.mark.xfail(reason="Implementation raises IndexError instead of input validation error when coverage column is missing", strict=False)
def test_run_rejects_coverage_without_depth_column(tmp_path: Path) -> None:
    # Intended behavior: reject malformed coverage input with a clear ValueError.
    fasta = _write_fasta(tmp_path / "ref.fasta", [("ref", "ACGT")])
    coverage = tmp_path / "coverage.tsv"
    coverage.write_text("position\n1\n2\n3\n4\n", encoding="utf-8")
    output = tmp_path / "boc.tsv"

    with pytest.raises(ValueError, match="coverage|column|malformed"):
        Boc(input=fasta, output=output, samplename="sample-bad", coverage=coverage).run()