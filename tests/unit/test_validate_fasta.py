from pathlib import Path
from typing import Callable, Generator

import pytest

from ViroConstrictor.validatefasta import (
    ContainsAmbiguities,
    ContainsSpecials,
    IsValidFasta,
    IsValidRef,
)

DATA_PATH = Path("tests/unit/data/temp/")


@pytest.mark.parametrize(
    "sequence, expected",
    [
        # Valid sequences without special characters
        ("ACTG", False),
        ("acgt", False),
        ("ACTactg", False),
        ("ACUGN", False),
        ("ACTGactg-", False),
        ("MRWSYKVHDBN", False),
        ("mrwsykvhdbn", False),
        # Valid sequence with ambiguity codes and dashes
        ("ACTG-MRWSY-KVHDBN", False),
        # Invalid sequences with special characters
        ("ACTG!", True),
        ("ACTG_ATCG", True),
        ("AC#TG", True),
        ("ACTG+ATCG", True),
        ("AC TG", True),  # Space is a special character
        ("ACTG*", True),  # '*' is considered special
        ("ACTG\nGTAC", True),  # Newline is special
        # Edge cases
        ("", False),  # Empty string has no specials
        ("-", False),  # Just a dash is valid
        ("N", False),  # Just ambiguity code is valid
    ],
)
def test_contains_specials(sequence: str, expected: bool) -> None:
    """Test the ContainsSpecials function to ensure it correctly identifies sequences with special characters."""

    assert ContainsSpecials(sequence) is expected


@pytest.mark.parametrize(
    "sequence, expected",
    [
        # Test sequences without ambiguities
        ("ACTG", False),
        ("acgt", False),
        ("ATGC-", False),
        ("", False),
        # Test sequences with ambiguities - lowercase
        ("actgu", True),  # u is an ambiguity
        ("actgm", True),  # m is an ambiguity (A or C)
        ("actgr", True),  # r is an ambiguity (A or G)
        ("actgw", True),  # w is an ambiguity (A or T)
        ("actgs", True),  # s is an ambiguity (G or C)
        ("actgy", True),  # y is an ambiguity (C or T)
        ("actgk", True),  # k is an ambiguity (G or T)
        ("actgv", True),  # v is an ambiguity (A, C, or G)
        ("actgh", True),  # h is an ambiguity (A, C, or T)
        ("actgd", True),  # d is an ambiguity (A, G, or T)
        ("actgb", True),  # b is an ambiguity (C, G, or T)
        ("actgn", True),  # n is an ambiguity (any base)
        # Test sequences with ambiguities - uppercase
        ("ACTGU", True),
        ("ACTGM", True),
        ("ACTGR", True),
        ("ACTGW", True),
        ("ACTGS", True),
        ("ACTGY", True),
        ("ACTGK", True),
        ("ACTGV", True),
        ("ACTGH", True),
        ("ACTGD", True),
        ("ACTGB", True),
        ("ACTGN", True),
        # Test mixed case sequences with ambiguities
        ("ACTGn", True),
        ("actgN", True),
        # Test sequences with multiple ambiguities
        ("ACNGT", True),
        ("ANNTG", True),
        ("NNNNN", True),
        # Test sequences with ambiguities and dashes
        ("A-N-G-T", True),
        ("-N-", True),
    ],
)
def test_contains_ambiguities(sequence: str, expected: bool) -> None:
    """Test the ContainsAmbiguities function to ensure it correctly identifies sequences with ambiguity characters."""
    assert ContainsAmbiguities(sequence) is expected


@pytest.fixture(scope="module")
def create_file_func() -> Generator[Callable, None, None]:

    DATA_PATH.mkdir(parents=True, exist_ok=True)

    created_files = []

    def create_file(file_name: str, content: str) -> Path:
        file_path = DATA_PATH / file_name
        with open(file_path, "w", encoding="utf-8") as f:
            f.write(content)
        created_files.append(file_path)
        return file_path

    yield create_file
    for file in created_files:
        if file.exists():
            file.unlink()
    DATA_PATH.rmdir()


@pytest.fixture(scope="module")
def temp_fasta_files(create_file_func) -> Generator[dict[Path, bool], None, None]:
    """
    Create temporary FASTA files for testing.
    """

    # Happy paths
    happy_files = [
        ("valid.fasta", ">seq1\nACTGACTG\n>seq2\nGTCAGTCA"),
        ("multi_line.fasta", ">seq1\nACGT\nGTCAGTCA\n>seq2\nGTCAGTCA"),
        ("ambiguous.fasta", ">seq1\nACNGACTG\n>seq2\nGTCAGTCA"),
        ("empty.fasta", ""),
    ]

    unhappy_files = [
        ("special.fasta", ">seq1\nACTG!ACTG\n>seq2\nGTCAGTCA"),
        ("mixed.fasta", ">seq1\nACTGACTG\n>seq2\nGTCAGTCA*\n>seq3\nACTGACTG"),
        (
            "missing_header.fasta",
            "ACTGACTG\n>seq2\nGTCAGTCA",
        ),
        (
            "inconsistent_line_breaks.fasta",
            ">seq1\nACTG\r>seq2\nGTCAGTCA\n>seq3\nACTGACTG\n>seq4\n",
        ),
    ]
    for file_name, content in happy_files + unhappy_files:
        create_file_func(file_name, content)

    happy_dict = {DATA_PATH / path: True for path, _ in happy_files}
    unhappy_dict = {DATA_PATH / path: False for path, _ in unhappy_files}
    fasta_dict = {**happy_dict, **unhappy_dict}
    yield fasta_dict


@pytest.fixture(scope="module")
def temp_ref_files(create_file_func) -> Generator[dict[Path, bool], None, None]:
    """Create temporary reference files for testing."""
    # Happy paths
    happy_files = [
        ("valid.fasta", ">seq1\nACTGACTG\n>seq2\nGTCAGTCA"),
        ("multi_line.fasta", ">seq1\nACGT\nGTCAGTCA\n>seq2\nGTCAGTCA"),
    ]

    unhappy_files = [
        ("empty.fasta", ""),
        ("ambiguous.fasta", ">seq1\nACNGACTG\n>seq2\nGTCAGTCA"),
        ("special.fasta", ">seq1\nACTG!ACTG\n>seq2\nGTCAGTCA"),
        ("mixed.fasta", ">seq1\nACTGACTG\n>seq2\nGTCAGTCA@\n>seq3\nACTGACTG"),
        (
            "missing_header.fasta",
            "ACTGACTG\n>seq2\nGTCAGTCA",
        ),
        (
            "inconsistent_line_breaks.fasta",
            ">seq1\nACTG\r>seq2\nGTCAGTCA\n>seq3\nACTGACTG\n>seq4\n",
        ),
    ]
    for file_name, content in happy_files + unhappy_files:
        create_file_func(file_name, content)

    happy_dict = {DATA_PATH / path: True for path, _ in happy_files}
    unhappy_dict = {DATA_PATH / path: False for path, _ in unhappy_files}
    fasta_dict = {**happy_dict, **unhappy_dict}
    yield fasta_dict


def test_is_valid_ref(temp_ref_files: dict[Path, bool]) -> None:
    """Test the IsValidRef function to ensure it correctly identifies valid reference files."""

    for file, is_valid in temp_ref_files.items():
        str_file = str(file)
        if is_valid:
            assert IsValidRef(str_file) is True
        else:
            assert IsValidRef(str_file) is False


def test_check_reference_file(temp_ref_files: dict[Path, bool]) -> None:
    """Test the IsValidRef function to ensure it correctly identifies valid reference files."""

    for file, is_valid in temp_ref_files.items():
        str_file = str(file)
        if is_valid:
            assert IsValidRef(str_file) is True
        else:
            assert IsValidRef(str_file) is False
