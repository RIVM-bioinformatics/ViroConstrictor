from pathlib import Path

import pytest

from ViroConstrictor.genbank import GenBank


def get_happy_cases() -> list[dict[str, Path | str]]:
    """Returns a list of valid test cases."""
    return [
        {
            "id": "test_normal",
            "input": Path("tests/unit/data/test_reference.gb"),
            "expected_fasta": Path("tests/unit/data/test_reference.fasta"),
            "expected_gff": Path("tests/unit/data/test_reference.gff"),
            "expected_target": "Influenza_A_virus",
        },
        {
            "id": "test_double",
            "input": Path("tests/unit/data/test_reference_double.gb"),
            "expected_fasta": Path("tests/unit/data/test_reference_double.fasta"),
            "expected_gff": Path("tests/unit/data/test_reference_double.gff"),
            "expected_target": "Influenza_A_virus",
        },
    ]


def get_unhappy_cases() -> list[dict[str, Path | str]]:
    """Returns a list of invalid test cases."""
    return [
        {
            "id": "test_bad_gb_file",
            "input": Path("tests/unit/data/test_reference_bad.gb"),
        }
    ]


def generate_test_cases() -> list[pytest.param]:
    happy_cases = [
        pytest.param(
            case["input"],
            case["expected_fasta"],
            case["expected_gff"],
            case["expected_target"],
            False,
            id=case["id"],
        )
        for case in get_happy_cases()
    ]
    unhappy_cases = [
        pytest.param(
            case["input"],
            None,
            None,
            None,
            True,
            id=case["id"],
        )
        for case in get_unhappy_cases()
    ]
    return happy_cases + unhappy_cases


@pytest.mark.parametrize(
    "path, expected_fasta, expected_gff, expected_target, should_raise",
    generate_test_cases(),
)
def test_split_genbank(
    path: Path,
    expected_fasta: Path | None,
    expected_gff: Path | None,
    expected_target: str | None,
    should_raise: bool,
):
    """Tests the split_genbank function with valid and invalid inputs."""
    if should_raise:
        with pytest.raises(ValueError, match="Error opening GenBank file:"):
            GenBank.split_genbank(path)
    else:
        actual_fasta, actual_gff, actual_target = GenBank.split_genbank(path)
        assert actual_fasta == expected_fasta
        assert actual_gff == expected_gff
        assert actual_target == expected_target
