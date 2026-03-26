"""
Unit tests for AggregateCombinedFiles script.

Tests the aggregation of combined tabular result files (mutations, coverage,
amplicon coverage) across multiple samples, including input validation, file
handling, and output format correctness.
"""

import sys
from argparse import ArgumentParser
from pathlib import Path

import pandas as pd
import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(PROJECT_ROOT.joinpath("ViroConstrictor/workflow")))
from ViroConstrictor.workflow.main.scripts.aggregate_combined_files import AggregateCombinedFiles  # isort:skip


def test_add_arguments_parses_valid_values() -> None:
    """The CLI parser should accept supported file types and default separator."""
    parser = ArgumentParser()
    AggregateCombinedFiles.add_arguments(parser)

    args = parser.parse_args(["--input", "in.txt", "--output", "out.txt", "--input_files", "a.tsv", "b.tsv", "--file_type", "mutations"])

    assert args.input == "in.txt"
    assert args.output == "out.txt"
    assert args.input_files == ["a.tsv", "b.tsv"]
    assert args.file_type == "mutations"
    assert args.separator == "\t"


def test_add_arguments_rejects_unknown_file_type() -> None:
    """The CLI parser should fail fast for unsupported aggregation types."""
    parser = ArgumentParser()
    AggregateCombinedFiles.add_arguments(parser)

    with pytest.raises(SystemExit):
        parser.parse_args(["--input", "in.txt", "--output", "out.txt", "--input_files", "a.tsv", "--file_type", "unknown_type"])


def test_run_delegates_to_aggregate_method(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """run() should trigger aggregation exactly once."""
    output = tmp_path / "aggregated.tsv"
    aggregator = AggregateCombinedFiles(input="unused", output=output, input_files=[], file_type="mutations")

    calls: list[str] = []

    def fake_aggregate() -> None:
        """Mock aggregation function for testing."""
        calls.append("called")

    monkeypatch.setattr(aggregator, "_aggregate_files", fake_aggregate)
    aggregator.run()

    assert calls == ["called"]


def test_run_aggregates_mutation_files_and_ignores_missing_or_empty_inputs(tmp_path: Path) -> None:
    """Only valid non-empty files should contribute rows to the output."""
    valid_a = tmp_path / "mutations_a.tsv"
    valid_b = tmp_path / "mutations_b.tsv"
    empty_file = tmp_path / "empty.tsv"
    missing_file = tmp_path / "does_not_exist.tsv"
    output = tmp_path / "all_mutations.tsv"

    valid_a.write_text("Sample\tVirus\tReference_ID\tPosition\tReference_Base\tVariant_Base\tDepth\nsample1\tVirusA\tRefA\t100\tA\tT\t50\n")
    valid_b.write_text("Sample\tVirus\tReference_ID\tPosition\tReference_Base\tVariant_Base\tDepth\nsample2\tVirusB\tRefB\t200\tC\tG\tNA\n")
    empty_file.write_text("")

    AggregateCombinedFiles(
        input="unused",
        output=output,
        input_files=[valid_a, missing_file, empty_file, valid_b],
        file_type="mutations",
    ).run()

    assert output.exists()
    result = pd.read_csv(output, sep="\t", dtype=str, keep_default_na=False)
    assert result.shape[0] == 2
    assert set(result["Sample"]) == {"sample1", "sample2"}
    assert result.loc[result["Sample"] == "sample2", "Depth"].item() == "NA"


def test_run_aggregates_amplicon_coverage_as_csv(tmp_path: Path) -> None:
    """Amplicon coverage inputs should be treated as comma-separated data."""
    amplicon_a = tmp_path / "amplicon_a.csv"
    amplicon_b = tmp_path / "amplicon_b.csv"
    output = tmp_path / "all_amplicon.csv"

    amplicon_a.write_text("Virus,Reference_ID,amplicon,mean_cov\nVirusA,RefA,amp1,100\n")
    amplicon_b.write_text("Virus,Reference_ID,amplicon,mean_cov\nVirusB,RefB,amp2,150\n")

    AggregateCombinedFiles(
        input="unused",
        output=output,
        input_files=[amplicon_a, amplicon_b],
        file_type="amplicon_coverage",
    ).run()

    assert output.exists()
    result = pd.read_csv(output)
    assert result.shape[0] == 2
    assert set(result["amplicon"]) == {"amp1", "amp2"}


def test_run_propagates_unexpected_parser_errors(tmp_path: Path) -> None:
    """Corrupt non-empty input should surface parsing failures instead of silent truncation."""
    corrupt = tmp_path / "corrupt.tsv"
    output = tmp_path / "all_mutations.tsv"

    # Unclosed quote forces a parser error for the default C engine.
    corrupt.write_text('Sample\tVirus\n"unterminated\tVirusA\n')

    aggregator = AggregateCombinedFiles(
        input="unused",
        output=output,
        input_files=[corrupt],
        file_type="mutations",
    )

    with pytest.raises(pd.errors.ParserError):
        aggregator.run()


@pytest.mark.parametrize(
    ("file_type", "expected"),
    [
        ("mutations", "Sample\tVirus\tReference_ID\tPosition\tReference_Base\tVariant_Base\tDepth\n"),
        (
            "coverage",
            "Sample_name\tVirus\tReference_ID\tWidth_at_mincov_1\tWidth_at_mincov_5\tWidth_at_mincov_10\tWidth_at_mincov_50\tWidth_at_mincov_100\n",
        ),
        ("amplicon_coverage", ""),
    ],
)
def test_run_writes_expected_empty_output_when_no_valid_input(tmp_path: Path, file_type: str, expected: str) -> None:
    """When no valid records exist, the output should still be initialized deterministically."""
    output = tmp_path / f"empty_{file_type}.txt"
    empty_input = tmp_path / "empty_input.txt"
    empty_input.write_text("")

    AggregateCombinedFiles(
        input="unused",
        output=output,
        input_files=[empty_input],
        file_type=file_type,
    ).run()

    assert output.read_text() == expected


def test_read_files_skips_empty_data_error(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Files that trigger EmptyDataError should be ignored while valid ones are kept."""
    broken = tmp_path / "broken.tsv"
    valid = tmp_path / "valid.tsv"
    broken.write_text("not_used")
    valid.write_text("also_not_used")

    aggregator = AggregateCombinedFiles(
        input="unused",
        output=tmp_path / "out.tsv",
        input_files=[],
        file_type="mutations",
    )

    def fake_read_csv(path: Path, sep: str = "\t", keep_default_na: bool = False, na_filter: bool = False) -> pd.DataFrame:
        """Mock CSV reader for testing edge cases."""
        if path == broken:
            raise pd.errors.EmptyDataError("empty")
        return pd.DataFrame([{"Sample": "sample1", "Depth": "9"}])

    monkeypatch.setattr(pd, "read_csv", fake_read_csv)
    frames = aggregator._read_files([broken, valid])

    assert len(frames) == 1
    assert frames[0].iloc[0]["Sample"] == "sample1"


@pytest.mark.xfail(reason="Current implementation does not validate schema for mutations/coverage files.", strict=True)
def test_intended_behavior_rejects_wrong_schema_for_mutation_aggregation(tmp_path: Path) -> None:
    """Intended behavior: malformed combined mutation files should be rejected."""
    malformed = tmp_path / "malformed_mutations.tsv"
    output = tmp_path / "all_mutations.tsv"

    malformed.write_text("foo\tbar\n1\t2\n")

    with pytest.raises(ValueError):
        AggregateCombinedFiles(
            input="unused",
            output=output,
            input_files=[malformed],
            file_type="mutations",
        ).run()


@pytest.mark.xfail(reason="Current implementation silently accepts unsupported file_type values.", strict=True)
def test_intended_behavior_rejects_unsupported_file_type_in_class_usage(tmp_path: Path) -> None:
    """Intended behavior: unsupported file types should raise a clear error."""
    output = tmp_path / "unsupported.txt"

    with pytest.raises(ValueError):
        AggregateCombinedFiles(
            input="unused",
            output=output,
            input_files=[],
            file_type="not_supported",
        ).run()
