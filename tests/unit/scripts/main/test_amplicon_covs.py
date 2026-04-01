"""
Unit tests for AmpliconCovs script.

Tests primer name parsing (supporting various formats and alt-allele notation),
BED primer file handling, coverage calculation per amplicon, and end-to-end
amplicon coverage aggregation.
"""

import sys
from argparse import ArgumentParser
from pathlib import Path

import pandas as pd
import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(PROJECT_ROOT.joinpath("ViroConstrictor/workflow")))
from ViroConstrictor.workflow.main.scripts.amplicon_covs import (  # noqa: E402, isort:skip
    AltName,
    AmpliconCovs,
    PrimerNameParser,
    ReadDirection,
)


def _write_bed(tmp_path: Path, lines: list[str]) -> Path:
    """
    Write BED format primer file to temporary path.

    Parameters
    ----------
    tmp_path : Path
        Temporary directory path.
    lines : list[str]
        List of BED format lines to write.

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
    Write per-position coverage file to temporary path.

    Parameters
    ----------
    tmp_path : Path
        Temporary directory path.
    values : list[int]
        List of coverage values per position (1-indexed).

    Returns
    -------
    Path
        Path to the written coverage file.
    """
    cov_path = tmp_path / "coverage.tsv"
    lines = [f"{idx + 1}\t{value}\n" for idx, value in enumerate(values)]
    cov_path.write_text("".join(lines), encoding="utf-8")
    return cov_path


@pytest.mark.parametrize(
    "token, expected",
    [
        ("alt", True),
        ("alt1", True),
        ("alternative", True),
        ("ALTLEFT", True),
        ("foo", False),
    ],
)
def test_alt_name_detection(token: str, expected: bool) -> None:
    """
    Verify that alt-allele notation is correctly detected in primer name tokens.

    Parametrized test that checks whether the AltName.is_valid_alt_name() method
    correctly identifies tokens that represent alternative alleles (e.g., "alt",
    "alt1", "alternative").
    """
    assert AltName.is_valid_alt_name(token) is expected


@pytest.mark.parametrize(
    "token, expected",
    [
        ("left", ReadDirection.FORWARD),
        ("FW", ReadDirection.FORWARD),
        ("reverse", ReadDirection.REVERSE),
        ("minus", ReadDirection.REVERSE),
    ],
)
def test_read_direction_from_string(token: str, expected: ReadDirection) -> None:
    """
    Verify that read direction strings are correctly converted to ReadDirection enum values.

    Parametrized test that validates the ReadDirection.from_string() method
    recognizes various string representations of forward and reverse directions.
    """
    assert ReadDirection.from_string(token) == expected


def test_read_direction_rejects_invalid_value() -> None:
    """
    Verify that invalid read direction strings are rejected with appropriate error handling.

    Tests both the is_valid_direction() validation method and the from_string() method
    to ensure invalid direction tokens raise ValueError with a descriptive message.
    """
    assert ReadDirection.is_valid_direction("sideways") is False
    with pytest.raises(ValueError, match="Unrecognized read direction"):
        ReadDirection.from_string("sideways")


@pytest.mark.parametrize(
    "primer_name, expected_name, expected_count, expected_alt, expected_direction",
    [
        (">ncov-2019_1_LEFT", "ncov-2019", 1, False, ReadDirection.FORWARD),
        ("ncov-2019_27_RIGHT_alt1", "ncov-2019", 27, True, ReadDirection.REVERSE),
        ("HAV_1_alt-IB_LEFT", "HAV", 1, True, ReadDirection.FORWARD),
        ("virus_strain_2_alt_RIGHT", "virus_strain", 2, True, ReadDirection.REVERSE),
        # Intended behavior: count should be amplicon number, not insert size.
        ("SARS-CoV-2_1200_100_LEFT_2", "SARS-CoV-2", 100, False, ReadDirection.FORWARD),
    ],
)
def test_parser_supported_formats(
    primer_name: str,
    expected_name: str,
    expected_count: int,
    expected_alt: bool,
    expected_direction: ReadDirection,
) -> None:
    """
    Verify that PrimerNameParser handles various supported primer name formats correctly.

    Parametrized test that validates the parser correctly extracts virus name, amplicon
    count, alt-allele status, and read direction from diverse primer naming conventions.
    """
    parsed = PrimerNameParser().parse(primer_name)
    assert parsed.name == expected_name
    assert parsed.count == expected_count
    assert parsed.alt is expected_alt
    assert parsed.direction == expected_direction


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
        ("SARS-COV-2_1_LEFT_1", True),
        ("SARS-COV-2_400_1_RIGHT_0", True),
        ("SARS-COV-2_1200_100_LEFT_2", True),
        # Invalid primer names
        ("wrong_name", False),  # Missing required components
        ("ncov-2019_LEFT", False),  # Missing primer number
        ("ncov-2019_alt", False),  # Missing direction
        ("ncov-2019_2_RIHTT", False),  # Misspelled direction
    ]

    error_keywords = ["Primer", "Unrecognized read direction", "does not match expected formats"]

    for primer_name, should_succeed in test_cases:
        if should_succeed:
            parsed = parser.parse(primer_name)
            assert isinstance(parsed.count, int)
            assert isinstance(parsed.direction, ReadDirection)
            assert parsed.original_string.lstrip(">") in primer_name
        else:
            with pytest.raises(ValueError) as exc_info:
                parser.parse(primer_name)
            assert any(keyword in str(exc_info.value) for keyword in error_keywords)


@pytest.mark.parametrize(
    "primer_name",
    [
        "wrong_name",
        "ncov-2019_LEFT",
        "ncov-2019_alt",
        "ncov-2019_2_RIHTT",
    ],
)
def test_parser_rejects_invalid_formats(primer_name: str) -> None:
    """
    Verify that PrimerNameParser rejects invalid primer name formats.

    Parametrized test that checks invalid primer names raise ValueError with
    appropriate error messages indicating parsing failures.
    """
    with pytest.raises(ValueError, match="Unrecognized read direction|does not match expected formats"):
        PrimerNameParser().parse(primer_name)


@pytest.mark.xfail(reason="Fallback path returns alt=None instead of bool False")
def test_fallback_parser_returns_boolean_alt_for_plain_three_part_name() -> None:
    """
    Test fallback parser return type consistency (XFAIL - intended behavior).

    Intended behavior: PrimerInfo.alt should always be a boolean value, but the
    fallback path currently returns None when parsing plain three-part names.
    """
    # Intended behavior: PrimerInfo.alt should always be boolean.
    parsed = PrimerNameParser()._fallback_parse("virus_8_LEFT")
    assert parsed.alt is False


def test_fallback_parser_handles_alt_in_either_position() -> None:
    """
    Test that fallback parser handles alt notation in different positions.

    Validates that the fallback parser correctly extracts name, count, and direction
    when alt notation appears as either the third or fourth component.
    """
    parser = PrimerNameParser()

    parsed_part3 = parser._fallback_parse("virus_3_alt1_RIGHT")
    assert parsed_part3.name == "virus"
    assert parsed_part3.count == 3
    assert parsed_part3.direction == ReadDirection.REVERSE

    parsed_part4 = parser._fallback_parse("virus_4_LEFT_alt2")
    assert parsed_part4.name == "virus"
    assert parsed_part4.count == 4
    assert parsed_part4.direction == ReadDirection.FORWARD


def test_fallback_parser_rejects_invalid_four_part_and_long_names() -> None:
    """
    Verify that fallback parser rejects invalid four-part and long primer names.

    Tests that the fallback parser raises ValueError when encountering malformed
    four-part names or names with more than four underscore-separated components.
    """
    parser = PrimerNameParser()

    with pytest.raises(ValueError, match="does not match expected formats"):
        parser._fallback_parse("virus_4_LEFT_notalt")

    with pytest.raises(ValueError, match="does not match expected formats"):
        parser._fallback_parse("virus_1_LEFT_alt_extra")


def test_add_arguments_registers_required_flags() -> None:
    """
    Verify that AmpliconCovs CLI parser registers all required command-line flags.

    Tests that the argument parser accepts input, output, coverages, and key flags
    and correctly stores their values when provided.
    """
    parser = ArgumentParser()
    AmpliconCovs.add_arguments(parser)
    args = parser.parse_args(["--input", "in.bed", "--output", "out.csv", "--coverages", "cov.tsv", "--key", "sampleX"])

    assert args.input == "in.bed"
    assert args.output == "out.csv"
    assert args.coverages == "cov.tsv"
    assert args.key == "sampleX"


def test_open_tsv_file_empty_returns_empty_dataframe(tmp_path: Path) -> None:
    """
    Verify that empty TSV files return an empty DataFrame.

    Tests that reading an empty file produces an empty DataFrame without errors.
    """
    empty_file = tmp_path / "empty.tsv"
    empty_file.write_text("", encoding="utf-8")
    assert AmpliconCovs._open_tsv_file(empty_file).empty


def test_open_tsv_file_raises_on_nan(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """
    Verify that _open_tsv_file raises ValueError when coverage contains NaN values.

    Tests that reading coverage files with NaN values triggers an error to prevent
    invalid coverage calculations.
    """
    bad_file = tmp_path / "bad.tsv"
    bad_file.write_text("1\t10\n", encoding="utf-8")

    def _fake_read_csv(*_args: object, **_kwargs: object) -> pd.DataFrame:
        return pd.DataFrame([[1, float("nan")]])

    monkeypatch.setattr(pd, "read_csv", _fake_read_csv)
    with pytest.raises(ValueError, match="contains NaN values"):
        AmpliconCovs._open_tsv_file(bad_file)


def test_calculate_amplicon_start_end_non_overlapping_boundaries(tmp_path: Path) -> None:
    """
    Verify correct amplicon boundary calculation for non-overlapping primer regions.

    Tests that start positions use the end of the LEFT primer and end positions
    use the start of the RIGHT primer for properly ordered primers.
    """
    bed_path = _write_bed(
        tmp_path,
        [
            "NC_045512.2\t10\t20\tvirus_1_LEFT\t0\t+\n",
            "NC_045512.2\t80\t90\tvirus_1_RIGHT\t0\t-\n",
            "NC_045512.2\t100\t110\tvirus_2_LEFT\t0\t+\n",
            "NC_045512.2\t180\t190\tvirus_2_RIGHT\t0\t-\n",
            "NC_045512.2\t200\t210\tvirus_3_LEFT\t0\t+\n",
            "NC_045512.2\t280\t290\tvirus_3_RIGHT\t0\t-\n",
        ],
    )
    covs = AmpliconCovs(input=bed_path, coverages=bed_path, key="sample", output=bed_path)
    primers = covs._split_primer_names(AmpliconCovs._open_tsv_file(bed_path))

    amplicons = covs._calculate_amplicon_start_end(primers)

    assert list(amplicons["amplicon_number"]) == [1, 2, 3]
    assert list(amplicons["start"]) == [20, 90, 190]
    assert list(amplicons["end"]) == [100, 200, 280]


def test_calculate_amplicon_start_end_falls_back_when_previous_reverse_missing(tmp_path: Path) -> None:
    """
    Verify fallback logic when previous amplicon lacks RIGHT primer.

    Tests that when a RIGHT primer is missing, the next amplicon's start position
    falls back to using the current amplicon's LEFT primer end position.
    """
    # Amplicon 1 has no RIGHT primer, so amplicon 2 start should use own LEFT end.
    bed_path = _write_bed(
        tmp_path,
        [
            "NC_045512.2\t10\t20\tvirus_1_LEFT\t0\t+\n",
            "NC_045512.2\t100\t110\tvirus_2_LEFT\t0\t+\n",
            "NC_045512.2\t170\t180\tvirus_2_RIGHT\t0\t-\n",
        ],
    )
    covs = AmpliconCovs(input=bed_path, coverages=bed_path, key="sample", output=bed_path)
    primers = covs._split_primer_names(AmpliconCovs._open_tsv_file(bed_path))

    amplicons = covs._calculate_amplicon_start_end(primers)
    second = amplicons.loc[amplicons["amplicon_number"] == 2].iloc[0]
    assert second["start"] == 110


def test_calculate_amplicon_start_end_falls_back_when_next_forward_missing(tmp_path: Path) -> None:
    """
    Verify fallback logic when next amplicon lacks LEFT primer.

    Tests that when a LEFT primer is missing, the current amplicon's end position
    falls back to using its own RIGHT primer start position.
    """
    # Amplicon 2 has no LEFT primer, so amplicon 1 end should use own RIGHT start.
    bed_path = _write_bed(
        tmp_path,
        [
            "NC_045512.2\t10\t20\tvirus_1_LEFT\t0\t+\n",
            "NC_045512.2\t80\t90\tvirus_1_RIGHT\t0\t-\n",
            "NC_045512.2\t170\t180\tvirus_2_RIGHT\t0\t-\n",
        ],
    )
    covs = AmpliconCovs(input=bed_path, coverages=bed_path, key="sample", output=bed_path)
    primers = covs._split_primer_names(AmpliconCovs._open_tsv_file(bed_path))

    amplicons = covs._calculate_amplicon_start_end(primers)
    first = amplicons.loc[amplicons["amplicon_number"] == 1].iloc[0]
    assert first["end"] == 80


def test_calculate_amplicon_start_end_raises_for_invalid_interval(tmp_path: Path) -> None:
    """
    Verify that invalid amplicon intervals (start >= end) raise ValueError.

    Tests that overlapping or inverted primer regions that result in start >= end
    are detected and rejected with an appropriate error message.
    """
    bed_path = _write_bed(
        tmp_path,
        [
            "NC_045512.2\t10\t100\tvirus_1_LEFT\t0\t+\n",
            "NC_045512.2\t50\t60\tvirus_1_RIGHT\t0\t-\n",
        ],
    )
    covs = AmpliconCovs(input=bed_path, coverages=bed_path, key="sample", output=bed_path)
    primers = covs._split_primer_names(AmpliconCovs._open_tsv_file(bed_path))

    with pytest.raises(ValueError, match="Invalid amplicon"):
        covs._calculate_amplicon_start_end(primers)


def test_calculate_mean_coverage_uses_one_based_inclusive_window() -> None:
    """
    Verify mean coverage calculation uses 1-based inclusive intervals.

    Tests that coverage values are correctly averaged over a 1-based inclusive
    position window (positions 2-4 inclusive).
    """
    coverages = pd.DataFrame({1: [10, 20, 30, 40, 50]})
    interval = pd.Series({"start": 2, "end": 4})
    assert AmpliconCovs._calculate_mean_coverage(interval, coverages) == pytest.approx(30.0)


def test_create_amplicon_names_list_pads_and_sorts(tmp_path: Path) -> None:
    """
    Verify that amplicon names are padded with zeros and sorted correctly.

    Tests that amplicon numbers are formatted with zero-padding for consistent
    alphanumeric sorting (e.g., virus_001, virus_012).
    """
    bed_path = _write_bed(
        tmp_path,
        [
            "NC_045512.2\t10\t20\tvirus_12_LEFT\t0\t+\n",
            "NC_045512.2\t80\t90\tvirus_1_RIGHT\t0\t-\n",
        ],
    )
    covs = AmpliconCovs(input=bed_path, coverages=bed_path, key="sample", output=bed_path)
    primers = covs._split_primer_names(AmpliconCovs._open_tsv_file(bed_path))
    assert covs._create_amplicon_names_list(primers) == ["virus_001", "virus_012"]


@pytest.mark.xfail(reason="Current implementation assumes index label 0 exists")
def test_create_amplicon_names_list_works_with_non_zero_based_index(tmp_path: Path) -> None:
    """
    Test amplicon name creation with non-standard DataFrame indices (XFAIL - intended behavior).

    Intended behavior: Amplicon names should be created correctly even when the DataFrame
    index is not zero-based, but the current implementation assumes index label 0 exists.
    """
    bed_path = _write_bed(
        tmp_path,
        [
            "NC_045512.2\t10\t20\tvirus_1_LEFT\t0\t+\n",
            "NC_045512.2\t80\t90\tvirus_2_RIGHT\t0\t-\n",
        ],
    )
    covs = AmpliconCovs(input=bed_path, coverages=bed_path, key="sample", output=bed_path)
    primers = covs._split_primer_names(AmpliconCovs._open_tsv_file(bed_path))
    shifted_index = primers.set_index(pd.Index([10, 11]))
    assert covs._create_amplicon_names_list(shifted_index) == ["virus_001", "virus_002"]


def test_run_end_to_end_produces_expected_coverage_values(tmp_path: Path) -> None:
    """
    Verify end-to-end amplicon coverage calculation produces correct values.

    Tests the complete workflow from primer BED and coverage file to output CSV,
    validating that mean coverage values are correctly computed for each amplicon.
    """
    bed_path = _write_bed(
        tmp_path,
        [
            "NC_045512.2\t1\t3\tvirus_1_LEFT\t0\t+\n",
            "NC_045512.2\t7\t9\tvirus_1_RIGHT\t0\t-\n",
            "NC_045512.2\t10\t12\tvirus_2_LEFT\t0\t+\n",
            "NC_045512.2\t16\t18\tvirus_2_RIGHT\t0\t-\n",
        ],
    )
    cov_path = _write_coverage(tmp_path, [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180])
    output_file = tmp_path / "amplicon_covs.csv"

    AmpliconCovs(input=bed_path, coverages=cov_path, key="sampleA", output=output_file).run()

    result = pd.read_csv(output_file, index_col=0)
    assert list(result.columns) == ["virus_001", "virus_002"]
    assert result.loc["sampleA", "virus_001"] == pytest.approx(65.0)
    assert result.loc["sampleA", "virus_002"] == pytest.approx(125.0)


def test_run_with_empty_primers_writes_only_index_row(tmp_path: Path) -> None:
    """
    Verify that empty primer BED file produces output with only sample index.

    Tests that when the primer BED file is empty, the output CSV contains only
    the sample name row with no amplicon columns.
    """
    bed_path = tmp_path / "empty.bed"
    bed_path.write_text("", encoding="utf-8")
    cov_path = _write_coverage(tmp_path, [10, 20, 30])
    output_file = tmp_path / "output.csv"

    AmpliconCovs(input=bed_path, coverages=cov_path, key="sample1", output=output_file).run()

    result = pd.read_csv(output_file, index_col=0)
    assert list(result.index) == ["sample1"]
    assert result.shape == (1, 0)
