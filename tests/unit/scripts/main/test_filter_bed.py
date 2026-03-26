"""
Unit tests for FilterBedInput script.

Tests BED file filtering by reference ID, preserving file format and literal
"NA" values while handling empty matches and missing files.
"""

import sys
from argparse import ArgumentParser
from pathlib import Path

import pandas as pd
import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(PROJECT_ROOT.joinpath("ViroConstrictor/workflow")))
from ViroConstrictor.workflow.main.scripts.filter_bed_input import FilterBedInput  # isort:skip


def test_add_arguments_parses_reference_id() -> None:
    """
    Verify that FilterBedInput CLI parser registers reference_id argument.

    Tests that the argument parser accepts input, output, and reference_id
    arguments and correctly stores their values.
    """
    parser = ArgumentParser()
    FilterBedInput.add_arguments(parser)
    args = parser.parse_args(["--input", "in.bed", "--output", "out.bed", "--reference_id", "REF_1"])
    assert args.reference_id == "REF_1"


def test_filter_bed_input_filters_matching_reference_rows(tmp_path: Path) -> None:
    """
    Verify BED rows are filtered to match the specified reference ID.

    Tests that only BED records with the first column (chrom) matching the
    specified reference_id are included in the output.
    """
    input_path = tmp_path / "input.bed"
    output_path = tmp_path / "output.bed"

    input_path.write_text(
        "REF_A\t10\t20\tamp1\t0\t+\nREF_B\t30\t40\tamp2\t0\t-\nREF_A\t50\t60\tamp3\t0\t+\n",
        encoding="utf-8",
    )

    FilterBedInput(input=str(input_path), output=str(output_path), reference_id="REF_A").run()

    df = pd.read_csv(output_path, sep="\t", header=None)
    assert df.shape == (2, 6)
    assert set(df.iloc[:, 0]) == {"REF_A"}
    assert set(df.iloc[:, 3]) == {"amp1", "amp3"}


def test_filter_bed_input_when_no_matches_writes_empty_file(tmp_path: Path) -> None:
    """
    Verify that non-matching reference produces empty output file.

    Tests that when no BED records match the specified reference_id, an empty
    output file is created.
    """
    input_path = tmp_path / "input.bed"
    output_path = tmp_path / "output.bed"

    input_path.write_text("REF_B\t10\t20\tamp1\t0\t+\n", encoding="utf-8")

    FilterBedInput(input=str(input_path), output=str(output_path), reference_id="REF_A").run()

    assert output_path.exists()
    assert output_path.read_text(encoding="utf-8") == ""


def test_filter_bed_input_preserves_literal_na_values(tmp_path: Path) -> None:
    """
    Verify that literal "NA" values in BED fields are preserved.

    Tests that the string "NA" appearing in BED file fields is not converted
    to NaN or treated as a missing value marker.
    """
    input_path = tmp_path / "input.bed"
    output_path = tmp_path / "output.bed"

    input_path.write_text("REF_A\t10\t20\tNA\t0\t+\n", encoding="utf-8")

    FilterBedInput(input=str(input_path), output=str(output_path), reference_id="REF_A").run()

    df = pd.read_csv(output_path, sep="\t", header=None, keep_default_na=False)
    assert df.iloc[0, 3] == "NA"


def test_filter_bed_input_propagates_missing_file_error(tmp_path: Path) -> None:
    """
    Verify that missing input file raises FileNotFoundError.

    Tests that the run() method fails with FileNotFoundError when the input
    BED file does not exist.
    """
    input_path = tmp_path / "does_not_exist.bed"
    output_path = tmp_path / "output.bed"

    with pytest.raises(FileNotFoundError):
        FilterBedInput(input=str(input_path), output=str(output_path), reference_id="REF_A").run()


def test_filter_bed_input_asserts_on_non_string_paths(tmp_path: Path) -> None:
    """
    Verify that non-string path types trigger AssertionError.

    Tests that the run() method uses assertions to validate that input and output
    are string types, failing fast on type mismatches.
    """
    input_path = tmp_path / "input.bed"
    input_path.write_text("REF_A\t10\t20\tamp1\t0\t+\n", encoding="utf-8")

    script = FilterBedInput(input=str(input_path), output=str(tmp_path / "out.bed"), reference_id="REF_A")
    script.input = [str(input_path)]  # type: ignore[assignment]
    with pytest.raises(AssertionError, match="Input"):
        script.run()


@pytest.mark.xfail(
    reason="Current script uses AssertionError for invalid path types; intended behavior is explicit ValueError with actionable message", strict=False
)
def test_filter_bed_input_rejects_invalid_types_with_clear_error(tmp_path: Path) -> None:
    """
    Test invalid path type rejection (XFAIL - intended behavior).

    Intended behavior: Invalid path types should raise ValueError with a clear
    message, but the current implementation uses AssertionError instead.
    """
    input_path = tmp_path / "input.bed"
    input_path.write_text("REF_A\t10\t20\tamp1\t0\t+\n", encoding="utf-8")

    script = FilterBedInput(input=str(input_path), output=str(tmp_path / "out.bed"), reference_id="REF_A")
    script.output = 123  # type: ignore[assignment]
    with pytest.raises(ValueError, match="output|string|path"):
        script.run()
