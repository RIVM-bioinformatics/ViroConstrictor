from contextlib import contextmanager
from pathlib import Path
from typing import Iterator

import pytest

from ViroConstrictor.workflow.main.scripts.amplicon_covs import AmpliconCovs


@contextmanager
def temporarily_modify_file(file_path: Path, new_line: str) -> Iterator[None]:
    """Context manager that temporarily modifies the file in a file for testing purposes."""
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
    """Test basic functionality of AmpliconCovs with valid input."""
    input_file = "tests/unit/data/ESIB_EQA_2024_SARS1_01_primers.bed"
    output_file = tmp_path / "ESIB_EQA_2024_SARS1_01_primers_output.csv"
    coverages_file = "tests/unit/data/ESIB_EQA_2024_SARS1_01_coverage.tsv"
    key = "ESIB_EQA_2024_SARS1_01"

    a = AmpliconCovs(input=input_file, output=output_file, key=key, coverages=coverages_file)
    a.run()

    assert output_file.exists()


def test_primer_name_validation(tmp_path: Path) -> None:
    """Test validation of primer names - both valid and invalid formats."""
    input_file = Path("tests/unit/data/ESIB_EQA_2024_SARS1_01_primers.bed")
    output_file = tmp_path / "test_output.csv"
    coverages_file = "tests/unit/data/ESIB_EQA_2024_SARS1_01_coverage.tsv"
    key = "ESIB_EQA_2024_SARS1_01"

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
        # Invalid primer names
        ("wrong_name", False),  # Missing required components
        ("ncov-2019_LEFT", False),  # Missing primer number
        ("ncov-2019_alt", False),  # Missing direction
        ("ncov-2019_2_RIHTT", False),  # Misspelled direction
    ]

    for primer_name, should_succeed in test_cases:
        bed_line = f"NC_045512.2\t1\t100\t{primer_name}\t0\t+\n"

        with temporarily_modify_file(input_file, bed_line):
            amplicon_covs = AmpliconCovs(input=input_file, output=output_file, key=key, coverages=coverages_file)

            if should_succeed:
                amplicon_covs.run()
                assert output_file.exists(), f"Failed for valid primer name: {primer_name}"
                output_file.unlink()  # Clean up for next iteration
            else:
                with pytest.raises(ValueError) as exc_info:
                    amplicon_covs.run()

                error_msg = str(exc_info.value)
                assert any(
                    keyword in error_msg for keyword in ["Primer", "Unrecognized read direction"]
                ), f"Unexpected error message for {primer_name}: {error_msg}"
