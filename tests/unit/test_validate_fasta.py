"""Tests for :mod:`ViroConstrictor.validatefasta`.

This module covers FASTA validation helpers, ambiguous sequence detection,
reference-header validation, and reference-file checks.
"""

from pathlib import Path

import pytest

from ViroConstrictor import validatefasta


def write_fasta(tmp_path: Path, name: str, content: str) -> Path:
    """Write a FASTA payload to a temporary file.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Temporary directory provided by pytest.
    name : str
        File name to create under ``tmp_path``.
    content : str
        FASTA content to write.

    Returns
    -------
    pathlib.Path
        Path to the written FASTA file.
    """
    fasta_file = tmp_path / name
    fasta_file.write_text(content, encoding="utf-8")
    return fasta_file


@pytest.mark.parametrize(
    "sequence, expected",
    [
        ("ACTG", False),
        ("acgt", False),
        ("ACTG-MRWSY-KVHDBN", False),
        ("", False),
        ("ACTG!", True),
        ("ACTG*", True),
        ("AC TG", True),
    ],
)
def test_contains_specials(sequence: str, expected: bool) -> None:
    """Test detection of special characters in sequences.

    Verifies that ContainsSpecials returns False for valid DNA bases and ambiguity codes,
    and True for special characters like punctuation or spaces.

    Parameters
    ----------
    sequence : str
        DNA sequence to test.
    expected : bool
        Whether the sequence contains special characters.
    """
    assert validatefasta.ContainsSpecials(sequence) is expected


@pytest.mark.parametrize(
    "sequence, expected",
    [
        ("ACTG", False),
        ("ATGC-", False),
        ("", False),
        ("ACTGN", True),
        ("actgm", True),
        ("A-N-G-T", True),
    ],
)
def test_contains_ambiguities(sequence: str, expected: bool) -> None:
    """Test detection of ambiguous bases in sequences.

    Verifies that ContainsAmbiguities returns False for unambiguous DNA bases and gaps,
    and True when ambiguity codes (N, M, etc.) are present.

    Parameters
    ----------
    sequence : str
        DNA sequence to test.
    expected : bool
        Whether the sequence contains ambiguous bases.
    """
    assert validatefasta.ContainsAmbiguities(sequence) is expected


def test_is_valid_fasta_accepts_none_keyword() -> None:
    """Test that IsValidFasta accepts 'NONE' as a keyword value."""
    assert validatefasta.IsValidFasta("NONE") is True


def test_is_valid_fasta_happy_path(tmp_path: Path) -> None:
    """Test IsValidFasta accepts well-formed FASTA files.

    Verifies that properly formatted FASTA files with valid sequences pass validation.

    Parameters
    ----------
    tmp_path : Path
        Temporary directory provided by pytest.
    """
    fasta = write_fasta(tmp_path, "valid.fasta", ">seq1\nACTGACTG\n>seq2\nGTCAGTCA\n")
    assert validatefasta.IsValidFasta(str(fasta)) is True


def test_is_valid_fasta_rejects_special_characters(tmp_path: Path) -> None:
    """Test IsValidFasta rejects FASTA files with special characters.

    Verifies that sequences containing invalid special characters fail validation.

    Parameters
    ----------
    tmp_path : Path
        Temporary directory provided by pytest.
    """
    fasta = write_fasta(tmp_path, "invalid.fasta", ">seq1\nACTG!ACTG\n")
    assert validatefasta.IsValidFasta(str(fasta)) is False


@pytest.mark.parametrize(
    "content, expected",
    [
        (">seq1\nACTGACTG\n", True),
        (">seq1\nACNT\n", False),
        ("this is not fasta content\n", False),
    ],
)
def test_is_valid_ref_file_driven_behavior(tmp_path: Path, content: str, expected: bool) -> None:
    """Test IsValidRef validates reference FASTA files with strict requirements.

    Verifies IsValidRef behavior with valid reference, ambiguous bases, and invalid format.

    Parameters
    ----------
    tmp_path : Path
        Temporary directory provided by pytest.
    content : str
        FASTA file content to test.
    expected : bool
        Whether the content should pass reference validation.
    """
    fasta = write_fasta(tmp_path, "reference.fasta", content)
    assert validatefasta.IsValidRef(str(fasta)) is expected


def test_check_reference_file_happy_path_no_errors(tmp_path: Path) -> None:
    """Test CheckReferenceFile accepts valid reference without errors.

    Verifies that properly formatted reference files with valid sequences pass without warnings.

    Parameters
    ----------
    tmp_path : Path
        Temporary directory provided by pytest.
    """
    fasta = write_fasta(tmp_path, "ref.fasta", ">valid_header\nACTGACTG\n")
    validatefasta.CheckReferenceFile(str(fasta), warnings_as_errors=False)


def test_check_reference_file_warns_on_ambiguities(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Test CheckReferenceFile warns on isolated ambiguous bases.

    Verifies that short stretches of ambiguous bases trigger warnings but not errors.

    Parameters
    ----------
    tmp_path : Path
        Temporary directory provided by pytest.
    monkeypatch : pytest.MonkeyPatch
        Pytest fixture for mocking behavior.
    """
    warnings: list[str] = []
    errors: list[str] = []
    monkeypatch.setattr(validatefasta.log, "warning", lambda message: warnings.append(message))
    monkeypatch.setattr(validatefasta.log, "error", lambda message: errors.append(message))

    fasta = write_fasta(tmp_path, "ref_warn.fasta", ">header\nACTGNN\n")
    validatefasta.CheckReferenceFile(str(fasta), warnings_as_errors=False)

    assert warnings
    assert errors == []


def test_check_reference_file_exits_on_long_ambiguity_stretch(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Test CheckReferenceFile exits on long stretches of ambiguous bases.

    Verifies that extended runs of ambiguous bases cause the function to exit with error.

    Parameters
    ----------
    tmp_path : Path
        Temporary directory provided by pytest.
    monkeypatch : pytest.MonkeyPatch
        Pytest fixture for mocking behavior.

    Raises
    ------
    SystemExit
        When long ambiguity stretch is detected.
    """
    errors: list[str] = []
    monkeypatch.setattr(validatefasta.log, "error", lambda message: errors.append(str(message)))

    fasta = write_fasta(tmp_path, "ref_error.fasta", ">header\nACTGNNNNN\n")
    with pytest.raises(SystemExit):
        validatefasta.CheckReferenceFile(str(fasta), warnings_as_errors=False)

    assert errors


def test_check_reference_file_warnings_become_errors(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Test CheckReferenceFile treats warnings as errors when flag is set.

    Verifies that with warnings_as_errors=True, even isolated ambiguities cause exit.

    Parameters
    ----------
    tmp_path : Path
        Temporary directory provided by pytest.
    monkeypatch : pytest.MonkeyPatch
        Pytest fixture for mocking behavior.

    Raises
    ------
    SystemExit
        When warnings_as_errors=True and any ambiguities are found.
    """
    errors: list[str] = []
    monkeypatch.setattr(validatefasta.log, "error", lambda message: errors.append(str(message)))

    fasta = write_fasta(tmp_path, "ref_warn_as_error.fasta", ">header\nACTGNN\n")
    with pytest.raises(SystemExit):
        validatefasta.CheckReferenceFile(str(fasta), warnings_as_errors=True)

    assert errors


def test_check_ref_header_returns_header_for_valid_input() -> None:
    """Test check_ref_header returns valid header unchanged."""
    assert validatefasta.check_ref_header("good_header") == "good_header"


@pytest.mark.parametrize(
    "header",
    ["", "CON", "...", "bad/header", "bad:header", "bad*header"],
)
def test_check_ref_header_rejects_invalid_inputs(header: str) -> None:
    """Test check_ref_header rejects invalid header formats.

    Verifies that empty, reserved, or character-restricted headers are rejected.

    Parameters
    ----------
    header : str
        Header string to validate.

    Raises
    ------
    SystemExit
        When header is invalid.
    """
    with pytest.raises(SystemExit):
        validatefasta.check_ref_header(header)
