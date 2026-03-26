import sys
from argparse import ArgumentParser
from pathlib import Path

import pandas as pd
import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(PROJECT_ROOT.joinpath("ViroConstrictor/workflow")))
from ViroConstrictor.workflow.main.scripts.amplicon_covs import (  # isort:skip
    AltName,
    AmpliconCovs,
    PrimerNameParser,
    ReadDirection,
)


def _write_bed(tmp_path: Path, lines: list[str]) -> Path:
    bed_path = tmp_path / "primers.bed"
    bed_path.write_text("".join(lines), encoding="utf-8")
    return bed_path


def _write_coverage(tmp_path: Path, values: list[int]) -> Path:
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
    assert ReadDirection.from_string(token) == expected


def test_read_direction_rejects_invalid_value() -> None:
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
    with pytest.raises(ValueError, match="Unrecognized read direction|does not match expected formats"):
        PrimerNameParser().parse(primer_name)


@pytest.mark.xfail(reason="Fallback path returns alt=None instead of bool False")
def test_fallback_parser_returns_boolean_alt_for_plain_three_part_name() -> None:
    # Intended behavior: PrimerInfo.alt should always be boolean.
    parsed = PrimerNameParser()._fallback_parse("virus_8_LEFT")
    assert parsed.alt is False


def test_fallback_parser_handles_alt_in_either_position() -> None:
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
    parser = PrimerNameParser()

    with pytest.raises(ValueError, match="does not match expected formats"):
        parser._fallback_parse("virus_4_LEFT_notalt")

    with pytest.raises(ValueError, match="does not match expected formats"):
        parser._fallback_parse("virus_1_LEFT_alt_extra")


def test_add_arguments_registers_required_flags() -> None:
    parser = ArgumentParser()
    AmpliconCovs.add_arguments(parser)
    args = parser.parse_args(["--input", "in.bed", "--output", "out.csv", "--coverages", "cov.tsv", "--key", "sampleX"])

    assert args.input == "in.bed"
    assert args.output == "out.csv"
    assert args.coverages == "cov.tsv"
    assert args.key == "sampleX"


def test_open_tsv_file_empty_returns_empty_dataframe(tmp_path: Path) -> None:
    empty_file = tmp_path / "empty.tsv"
    empty_file.write_text("", encoding="utf-8")
    assert AmpliconCovs._open_tsv_file(empty_file).empty


def test_open_tsv_file_raises_on_nan(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    bad_file = tmp_path / "bad.tsv"
    bad_file.write_text("1\t10\n", encoding="utf-8")

    def _fake_read_csv(*_args: object, **_kwargs: object) -> pd.DataFrame:
        return pd.DataFrame([[1, float("nan")]])

    monkeypatch.setattr(pd, "read_csv", _fake_read_csv)
    with pytest.raises(ValueError, match="contains NaN values"):
        AmpliconCovs._open_tsv_file(bad_file)


def test_calculate_amplicon_start_end_non_overlapping_boundaries(tmp_path: Path) -> None:
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
    coverages = pd.DataFrame({1: [10, 20, 30, 40, 50]})
    interval = pd.Series({"start": 2, "end": 4})
    assert AmpliconCovs._calculate_mean_coverage(interval, coverages) == 30.0


def test_create_amplicon_names_list_pads_and_sorts(tmp_path: Path) -> None:
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
    assert result.loc["sampleA", "virus_001"] == 65.0
    assert result.loc["sampleA", "virus_002"] == 125.0


def test_run_with_empty_primers_writes_only_index_row(tmp_path: Path) -> None:
    bed_path = tmp_path / "empty.bed"
    bed_path.write_text("", encoding="utf-8")
    cov_path = _write_coverage(tmp_path, [10, 20, 30])
    output_file = tmp_path / "output.csv"

    AmpliconCovs(input=bed_path, coverages=cov_path, key="sample1", output=output_file).run()

    result = pd.read_csv(output_file, index_col=0)
    assert list(result.index) == ["sample1"]
    assert result.shape == (1, 0)
