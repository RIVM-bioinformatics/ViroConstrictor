import sys
from contextlib import contextmanager
from pathlib import Path
from typing import Iterator

import pytest
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PROJECT_ROOT.joinpath("ViroConstrictor/workflow")))
from ViroConstrictor.workflow.main.scripts.amplicon_covs import (  # isort:skip
    AmpliconCovs,
    PrimerNameParser,
    ReadDirection,
)


@contextmanager
def temporarily_modify_file(file_path: Path, new_line: str) -> Iterator[None]:
    """
    Temporarily replace a file's contents and restore them on exit.

    Parameters
    ----------
    file_path : Path
        Path to the file to modify.
    new_line : str
        New content to write to the file while inside the context.

    Yields
    ------
    None
        Control back to the context block. Original file contents are restored on exit.
    """
    with open(file_path, "r", encoding="utf-8") as file:
        lines = file.readlines()

    try:
        with open(file_path, "w", encoding="utf-8") as file:
            file.writelines(new_line)
        yield

    finally:
        # Restore original lines
        with open(file_path, "w", encoding="utf-8") as file:
            file.writelines(lines)


def test_amplicon_covs(tmp_path: Path) -> None:
    """
    Validate basic functionality of `AmpliconCovs` using sample BED and coverage files.

    Parameters
    ----------
    tmp_path : Path
        Temporary directory fixture provided by pytest for output files.

    Returns
    -------
    None
    """
    input_file = PROJECT_ROOT / "tests" / "unit" / "data" / "ESIB_EQA_2024_SARS1_01_primers.bed"
    output_file = tmp_path / "ESIB_EQA_2024_SARS1_01_primers_output.csv"
    coverages_file = PROJECT_ROOT / "tests" / "unit" / "data" / "ESIB_EQA_2024_SARS1_01_coverage.tsv"
    key = "ESIB_EQA_2024_SARS1_01"

    a = AmpliconCovs(input=input_file, output=output_file, key=key, coverages=coverages_file)
    a.run()

    assert output_file.exists()


def test_primer_name_validation() -> None:
    """
    Verify the `PrimerNameParser` accepts valid primer name formats and rejects invalid ones.

    The test iterates over several example primer names and checks that parsing
    either succeeds (returning an object with an integer `count`) or raises
    a `ValueError` for invalid formats.

    Returns
    -------
    None
    """
    parser = PrimerNameParser()

    # Test cases: (primer_name, should_succeed)
    test_cases = [
        # Valid primer names
        (">ncov-2019_1_LEFT", True),
        (">ncov-2019_27_RIGHT_alt1", True),
        (">HAV_1_alt-IB_LEFT", True),
        (">HAV_2_alt-IIIA_RIGHT", True),
        ("MBTuni_1_alt1_LEFT", True),
        ("MBTuni_1_RIGHT", True),
        ("MeV_2_RIGHT", True),
        ("MeV_2_alt2_RIGHT", True),
        (">MuV-NGS_1_alt21_LEFT", True),
        ("MuV-NGS_0_RIGHT", True),
        ("ncov-2019_2_right", True),
        ("ncov-2019_2_LEFt", True),
        ("ncov-2019_2_RIGHT", True),
        ("ncov-2019_2_alt_LEFT", True),
        ("ncov-2019_2_right_ALT1", True),
        ("virus_strain_1_LEFT", True),
        ("virus_strain_1_alt_RIGHT", True),
        ("MuV-NGS_19-alt22_LEFT", True),
        ("SARS-CoV-2_12_RIGHT_0", True),
        # Invalid primer names
        ("wrong_name", False),  # Missing required components
        ("ncov-2019_LEFT", False),  # Missing primer number
        ("ncov-2019_alt", False),  # Missing direction
        ("ncov-2019_2_RIHTT", False),  # Misspelled direction
    ]

    for primer_name, should_succeed in test_cases:
        if should_succeed:
            info = parser.parse(primer_name)
            assert info.original_string.lstrip(">") in primer_name
            assert isinstance(info.count, int)
        else:
            with pytest.raises(ValueError) as exc_info:
                parser.parse(primer_name)

            error_msg = str(exc_info.value)
            assert any(
                keyword in error_msg for keyword in ["Primer", "Unrecognized read direction", "does not match expected formats"]
            ), f"Unexpected error message for {primer_name}: {error_msg}"


def _write_bed(tmp_path: Path, lines: list[str]) -> Path:
    """
    Write primer BED lines to a file in the provided temporary directory.

    Parameters
    ----------
    tmp_path : Path
        Temporary directory to write the file into.
    lines : list[str]
        Lines to write to the BED file (each including trailing newline).

    Returns
    -------
    Path
        Path to the written BED file.
    """

    bed_path = tmp_path / "primers.bed"
    bed_path.write_text("".join(lines), encoding="utf-8")
    return bed_path


def _write_coverage(tmp_path: Path, values: list[int]) -> Path:
    """
    Create a simple coverage TSV file with one coverage value per line.

    Parameters
    ----------
    tmp_path : Path
        Temporary directory to write the coverage file into.
    values : list[int]
        Coverage values to write; each value is written on its own row with
        a 1-based index in the first column.

    Returns
    -------
    Path
        Path to the written coverage TSV file.
    """

    cov_path = tmp_path / "coverage.tsv"
    lines = [f"{idx + 1}\t{value}\n" for idx, value in enumerate(values)]
    cov_path.write_text("".join(lines), encoding="utf-8")
    return cov_path


def test_open_tsv_file_empty(tmp_path: Path) -> None:
    """
    Ensure `_open_tsv_file` returns an empty DataFrame for an empty file.

    Parameters
    ----------
    tmp_path : Path
        Pytest temporary directory fixture.

    Returns
    -------
    None
    """

    empty_file = tmp_path / "empty.tsv"
    empty_file.write_text("", encoding="utf-8")

    df = AmpliconCovs._open_tsv_file(empty_file)

    assert df.empty


def test_open_tsv_file_nan_raises(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """
    Verify that `_open_tsv_file` raises a `ValueError` when the file contains NaN.

    Parameters
    ----------
    monkeypatch : pytest.MonkeyPatch
        Pytest fixture for temporarily patching attributes.
    tmp_path : Path
        Temporary directory fixture for creating the bad TSV file.

    Returns
    -------
    None
    """

    bad_file = tmp_path / "bad.tsv"
    bad_file.write_text("1\t10\n", encoding="utf-8")
    bad_df = pd.DataFrame([[1, float("nan")]])

    def _fake_read_csv(*_args: object, **_kwargs: object) -> pd.DataFrame:
        return bad_df

    monkeypatch.setattr(pd, "read_csv", _fake_read_csv)

    with pytest.raises(ValueError, match="contains NaN values"):
        AmpliconCovs._open_tsv_file(bad_file)


def test_calculate_amplicon_start_end_two_amplicons(tmp_path: Path) -> None:
    """
    Calculate start/end coordinates for two amplicons from paired primer entries.

    Parameters
    ----------
    tmp_path : Path
        Temporary directory fixture for creating a temporary BED file.

    Returns
    -------
    None
    """

    bed_lines = [
        "NC_045512.2\t10\t20\tvirus_1_LEFT\t0\t+\n",
        "NC_045512.2\t80\t90\tvirus_1_RIGHT\t0\t-\n",
        "NC_045512.2\t100\t110\tvirus_2_LEFT\t0\t+\n",
        "NC_045512.2\t170\t180\tvirus_2_RIGHT\t0\t-\n",
    ]
    bed_path = _write_bed(tmp_path, bed_lines)
    primers = AmpliconCovs._open_tsv_file(bed_path)
    covs = AmpliconCovs(input=bed_path, coverages=bed_path, key="sample", output=bed_path)
    primers = covs._split_primer_names(primers)

    amplicon_sizes = covs._calculate_amplicon_start_end(primers)

    first = amplicon_sizes.loc[amplicon_sizes["amplicon_number"] == 1].iloc[0]
    second = amplicon_sizes.loc[amplicon_sizes["amplicon_number"] == 2].iloc[0]
    assert (first["start"], first["end"]) == (20, 100)
    assert (second["start"], second["end"]) == (90, 170)


def test_calculate_amplicon_start_end_invalid_raises(tmp_path: Path) -> None:
    """
    Ensure `_calculate_amplicon_start_end` raises for overlapping/invalid primer pairs.

    Parameters
    ----------
    tmp_path : Path
        Temporary directory fixture for creating a temporary BED file.

    Returns
    -------
    None
    """

    bed_lines = [
        "NC_045512.2\t10\t100\tvirus_1_LEFT\t0\t+\n",
        "NC_045512.2\t50\t60\tvirus_1_RIGHT\t0\t-\n",
    ]
    bed_path = _write_bed(tmp_path, bed_lines)
    primers = AmpliconCovs._open_tsv_file(bed_path)
    covs = AmpliconCovs(input=bed_path, coverages=bed_path, key="sample", output=bed_path)
    primers = covs._split_primer_names(primers)

    with pytest.raises(ValueError, match="Invalid amplicon"):
        covs._calculate_amplicon_start_end(primers)


def test_calculate_mean_coverage() -> None:
    """
    Validate mean coverage calculation over a given start/end interval.

    Returns
    -------
    None
    """

    coverages = AmpliconCovs._open_tsv_file(
        Path(__file__).resolve().parents[3] / "tests" / "unit" / "data" / "ESIB_EQA_2024_SARS1_01_coverage.tsv",
        index_col=0,
    )
    input_array = pd.Series({"start": 2, "end": 4})
    mean_cov = AmpliconCovs._calculate_mean_coverage(input_array, coverages)

    assert mean_cov == round(float(coverages.iloc[1:4].mean().values[0]), 2)


def test_create_amplicon_names_list(tmp_path: Path) -> None:
    """
    Create normalized amplicon name list from primer entries (zero-padded numbers).

    Parameters
    ----------
    tmp_path : Path
        Temporary directory fixture for creating the BED file used in the test.

    Returns
    -------
    None
    """

    bed_lines = [
        "NC_045512.2\t10\t20\tvirus_1_LEFT\t0\t+\n",
        "NC_045512.2\t80\t90\tvirus_12_RIGHT\t0\t-\n",
    ]
    bed_path = _write_bed(tmp_path, bed_lines)
    primers = AmpliconCovs._open_tsv_file(bed_path)
    covs = AmpliconCovs(input=bed_path, coverages=bed_path, key="sample", output=bed_path)
    primers = covs._split_primer_names(primers)

    amplicon_names = covs._create_amplicon_names_list(primers)

    assert amplicon_names == ["virus_001", "virus_012"]


def test_amplicon_covs_empty_primers(tmp_path: Path) -> None:
    """
    Run `AmpliconCovs` with an empty primer BED file to ensure output is produced.

    Parameters
    ----------
    tmp_path : Path
        Temporary directory fixture for creating input and output files.

    Returns
    -------
    None
    """

    bed_path = tmp_path / "empty.bed"
    bed_path.write_text("", encoding="utf-8")
    cov_path = _write_coverage(tmp_path, [10, 20, 30])
    output_file = tmp_path / "output.csv"

    covs = AmpliconCovs(input=bed_path, coverages=cov_path, key="sample1", output=output_file)
    covs.run()

    assert output_file.exists()
    assert "sample1" in output_file.read_text(encoding="utf-8")


def test_primer_name_parser_fallback_errors() -> None:
    """
    Confirm `PrimerNameParser` raises `ValueError` for inputs that don't match any
    known primer name formats.

    Returns
    -------
    None
    """

    parser = PrimerNameParser()

    with pytest.raises(ValueError, match="does not match expected formats"):
        parser.parse("bad")

    with pytest.raises(ValueError, match="does not match expected formats"):
        parser.parse("ncov-2019_1")


def test_read_direction_from_string() -> None:
    """
    Verify `ReadDirection.from_string` parses common direction labels and
    raises for unrecognized values.

    Returns
    -------
    None
    """

    assert ReadDirection.from_string("LEFT") == ReadDirection.FORWARD
    assert ReadDirection.from_string("reverse") == ReadDirection.REVERSE

    with pytest.raises(ValueError, match="Unrecognized read direction"):
        ReadDirection.from_string("sideways")
