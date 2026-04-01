"""
Unit tests for ConcatAmpliconCovs script.

Tests concatenation of per-sample amplicon coverage CSV files with proper
index alignment, missing value handling, and outer-join merging.
"""

import sys
from argparse import ArgumentParser
from pathlib import Path

import pandas as pd
import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(PROJECT_ROOT.joinpath("ViroConstrictor/workflow")))
from ViroConstrictor.workflow.main.scripts.concat_amplicon_covs import ConcatAmpliconCovs  # noqa: E402, isort:skip


def _write_amplicon_csv(path: Path, rows: list[tuple[str, float, float]]) -> Path:
    """
    Write amplicon coverage CSV file to path.

    Parameters
    ----------
    path : Path
        Destination file path.
    rows : list[tuple[str, float, float]]
        List of (amplicon_name, mean_coverage, median_coverage) tuples.

    Returns
    -------
    Path
        Path to the written CSV file.
    """
    df = pd.DataFrame(rows, columns=["amplicon_names", "mean_cov", "median_cov"])
    df.to_csv(path, index=False)
    return path


def test_add_arguments_parses_input_coverages_list() -> None:
    """
    Verify that ConcatAmpliconCovs CLI parser registers input_coverages list.

    Tests that the argument parser accepts input, output, and input_coverages
    list arguments and correctly stores their values.
    """
    parser = ArgumentParser()
    ConcatAmpliconCovs.add_arguments(parser)

    args = parser.parse_args(["--input", "unused", "--output", "combined.csv", "--input_coverages", "a.csv", "b.csv"])

    assert args.input_coverages == ["a.csv", "b.csv"]


def test_read_csv_uses_amplicon_names_as_index(tmp_path: Path) -> None:
    """
    Verify that amplicon names are used as DataFrame index.

    Tests that _read_csv() reads the amplicon_names column and sets it as the
    DataFrame index, with coverage columns as remaining columns.
    """
    csv_path = _write_amplicon_csv(tmp_path / "sample.csv", [("amp_1", 10.0, 9.0), ("amp_2", 11.0, 10.0)])

    df = ConcatAmpliconCovs._read_csv(csv_path)

    assert list(df.index) == ["amp_1", "amp_2"]
    assert list(df.columns) == ["mean_cov", "median_cov"]


def test_run_concatenates_multiple_files_and_fills_missing_values(tmp_path: Path) -> None:
    """
    Verify multiple files are concatenated with outer join and missing values filled.

    Tests that the run() method performs an outer join on multiple input coverage
    files and fills missing values with NaN where amplicons are not present in all files.
    """
    first = tmp_path / "first.csv"
    second = tmp_path / "second.csv"
    output = tmp_path / "combined.csv"

    pd.DataFrame(
        [
            {"amplicon_names": "amp_1", "mean_cov": 10.0, "median_cov": 9.0},
            {"amplicon_names": "amp_2", "mean_cov": 20.0, "median_cov": 19.0},
        ]
    ).to_csv(first, index=False)
    pd.DataFrame(
        [
            {"amplicon_names": "amp_3", "mean_cov": 30.0, "extra_stat": 3.3},
        ]
    ).to_csv(second, index=False)

    ConcatAmpliconCovs(input=[""], input_coverages=[first, second], output=output).run()

    merged = pd.read_csv(output, index_col=0)
    assert "mean_cov" in merged.columns
    assert "median_cov" in merged.columns
    assert "extra_stat" in merged.columns
    assert pd.isna(merged.loc["amp_3", "median_cov"])
    assert float(merged.loc["amp_1", "mean_cov"]) == 10.0


def test_run_supports_single_input_when_input_is_string(tmp_path: Path) -> None:
    """
    Verify single input handling when input is a string instead of list.

    Tests that the run() method correctly processes a single coverage file
    when the input attribute is set to a string path rather than a list.
    """
    input_csv = _write_amplicon_csv(tmp_path / "single.csv", [("amp_10", 100.0, 99.0)])
    output = tmp_path / "combined.csv"

    script = ConcatAmpliconCovs(input=[""], input_coverages=[input_csv], output=output)
    script.input = str(input_csv)
    script.run()

    merged = pd.read_csv(output, index_col=0)
    assert list(merged.index) == ["amp_10"]
    assert float(merged.loc["amp_10", "mean_cov"]) == 100.0


def test_run_raises_for_missing_input_file(tmp_path: Path) -> None:
    """
    Verify that missing input file raises FileNotFoundError.

    Tests that the run() method fails with FileNotFoundError when an input
    coverage file does not exist.
    """
    output = tmp_path / "combined.csv"

    with pytest.raises(FileNotFoundError):
        ConcatAmpliconCovs(input=[""], input_coverages=[tmp_path / "missing.csv"], output=output).run()


def test_read_csv_raises_when_amplicon_names_column_missing(tmp_path: Path) -> None:
    """
    Verify that missing amplicon_names column raises KeyError.

    Tests that _read_csv() raises KeyError when the required amplicon_names
    column is not present in the CSV file.
    """
    bad = tmp_path / "bad.csv"
    pd.DataFrame([{"name": "amp_1", "mean_cov": 10.0}]).to_csv(bad, index=False)

    with pytest.raises(KeyError):
        ConcatAmpliconCovs._read_csv(bad)


@pytest.mark.xfail(reason="Empty input list propagates pandas concat error instead of clear validation error", strict=False)
def test_run_rejects_empty_input_list_with_clear_error(tmp_path: Path) -> None:
    """
    Test empty input list rejection (XFAIL - intended behavior).

    Intended behavior: Empty input_coverages list should raise ValueError with
    an explicit message, but the current implementation propagates pandas concat error.
    """
    # Intended behavior: raise ValueError with an explicit message before calling pandas.concat.
    output = tmp_path / "combined.csv"

    with pytest.raises(ValueError, match="input|coverage|empty"):
        ConcatAmpliconCovs(input=[""], input_coverages=[], output=output).run()
