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
    assert validatefasta.ContainsAmbiguities(sequence) is expected


def test_is_valid_fasta_accepts_none_keyword() -> None:
    assert validatefasta.IsValidFasta("NONE") is True


def test_is_valid_fasta_happy_path(tmp_path: Path) -> None:
    fasta = write_fasta(tmp_path, "valid.fasta", ">seq1\nACTGACTG\n>seq2\nGTCAGTCA\n")
    assert validatefasta.IsValidFasta(str(fasta)) is True


def test_is_valid_fasta_rejects_special_characters(tmp_path: Path) -> None:
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
    fasta = write_fasta(tmp_path, "reference.fasta", content)
    assert validatefasta.IsValidRef(str(fasta)) is expected


def test_check_reference_file_happy_path_no_errors(tmp_path: Path) -> None:
    fasta = write_fasta(tmp_path, "ref.fasta", ">valid_header\nACTGACTG\n")
    validatefasta.CheckReferenceFile(str(fasta), warnings_as_errors=False)


def test_check_reference_file_warns_on_ambiguities(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    warnings: list[str] = []
    errors: list[str] = []
    monkeypatch.setattr(validatefasta.log, "warning", lambda message: warnings.append(message))
    monkeypatch.setattr(validatefasta.log, "error", lambda message: errors.append(message))

    fasta = write_fasta(tmp_path, "ref_warn.fasta", ">header\nACTGNN\n")
    validatefasta.CheckReferenceFile(str(fasta), warnings_as_errors=False)

    assert warnings
    assert errors == []


def test_check_reference_file_exits_on_long_ambiguity_stretch(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    errors: list[str] = []
    monkeypatch.setattr(validatefasta.log, "error", lambda message: errors.append(str(message)))

    fasta = write_fasta(tmp_path, "ref_error.fasta", ">header\nACTGNNNNN\n")
    with pytest.raises(SystemExit):
        validatefasta.CheckReferenceFile(str(fasta), warnings_as_errors=False)

    assert errors


def test_check_reference_file_warnings_become_errors(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    errors: list[str] = []
    monkeypatch.setattr(validatefasta.log, "error", lambda message: errors.append(str(message)))

    fasta = write_fasta(tmp_path, "ref_warn_as_error.fasta", ">header\nACTGNN\n")
    with pytest.raises(SystemExit):
        validatefasta.CheckReferenceFile(str(fasta), warnings_as_errors=True)

    assert errors


def test_check_ref_header_returns_header_for_valid_input() -> None:
    assert validatefasta.check_ref_header("good_header") == "good_header"


@pytest.mark.parametrize(
    "header",
    ["", "CON", "...", "bad/header", "bad:header", "bad*header"],
)
def test_check_ref_header_rejects_invalid_inputs(header: str) -> None:
    with pytest.raises(SystemExit):
        validatefasta.check_ref_header(header)
