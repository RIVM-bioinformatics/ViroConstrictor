"""Unit tests for GenBank parsing and split helpers."""

import shutil
from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from ViroConstrictor.genbank import GenBank

PROJECT_ROOT = Path(__file__).resolve().parents[2]
UNIT_DATA_DIR = PROJECT_ROOT / "tests" / "unit" / "data"


@pytest.mark.parametrize(
    "file_name, expected",
    [
        ("reference.gb", True),
        ("reference.gbk", True),
        ("reference.genbank", True),
        ("reference.fasta", False),
    ],
)
def test_is_genbank_by_extension(file_name: str, expected: bool) -> None:
    """Test GenBank file extension detection.

    Verifies that GenBank.is_genbank() correctly identifies GenBank file
    extensions (.gb, .gbk, .genbank) as valid, while rejecting other
    extensions like .fasta.

    Parameters
    ----------
    file_name : str
        Name of the file to test (parametrized).
    expected : bool
        Expected result from is_genbank() (parametrized).
    """
    assert GenBank.is_genbank(Path(file_name)) is expected


def test_open_genbank_reads_valid_file() -> None:
    """Test parsing a valid GenBank file.

    Verifies that GenBank.open_genbank() successfully reads a GenBank file
    and returns a list of SeqRecord objects containing organism annotations.
    """
    records = GenBank.open_genbank(UNIT_DATA_DIR / "test_reference.gb")

    assert len(records) == 1
    assert records[0].annotations.get("organism")


def test_open_genbank_rejects_non_genbank_extension(tmp_path: Path) -> None:
    """Test that open_genbank rejects non-GenBank files.

    Verifies that GenBank.open_genbank() raises a ValueError with descriptive
    error message ("is not a GenBank file") when given a file with an
    unsupported extension.

    Parameters
    ----------
    tmp_path : Path
        Temporary directory provided by pytest.
    """
    fasta_path = tmp_path / "not_genbank.fasta"
    fasta_path.write_text(">id\nACGT\n", encoding="utf-8")

    with pytest.raises(ValueError, match=r"is not a GenBank file"):
        GenBank.open_genbank(fasta_path)


def test_open_genbank_wraps_parse_errors() -> None:
    """Test that open_genbank wraps parser errors.

    Verifies that GenBank.open_genbank() catches Biopython parsing exceptions
    and re-raises them as ValueError with context about the file being
    parsed.
    """
    with pytest.raises(ValueError, match=r"Error opening GenBank file"):
        GenBank.open_genbank(UNIT_DATA_DIR / "test_reference_bad.gb")


def test_parse_target_normalizes_organism_name() -> None:
    """Test that _parse_target normalizes organism names.

    Verifies that GenBank._parse_target() extracts organism names from
    annotations, normalizes them by removing strain details (e.g., "H1N1"),
    and replaces spaces with underscores for safe filenames.
    """
    records = [
        SeqRecord(Seq("ACGT"), id="r1", annotations={"organism": "Influenza A virus (H1N1)"}),
        SeqRecord(Seq("ACGT"), id="r2", annotations={"organism": "Influenza A virus (H3N2)"}),
    ]

    assert GenBank._parse_target(records) == "Influenza_A_virus"


def test_parse_target_rejects_dissimilar_organisms() -> None:
    """Test that _parse_target rejects dissimilar organisms.

    Verifies that GenBank._parse_target() raises ValueError when organism
    annotations across records do not meet the required similarity threshold
    (must be similar enough to represent a single target).
    """
    records = [
        SeqRecord(Seq("ACGT"), id="r1", annotations={"organism": "Influenza A virus"}),
        SeqRecord(Seq("ACGT"), id="r2", annotations={"organism": "Human adenovirus"}),
    ]

    with pytest.raises(ValueError, match=r"sufficiently similar organism annotations"):
        GenBank._parse_target(records)


def test_parse_target_ignores_records_without_organism() -> None:
    """Test that _parse_target ignores records without organism annotations.

    Verifies that GenBank._parse_target() gracefully handles records missing
    organism annotations by ignoring them and deriving the target name from
    records that have complete organism information.
    """
    records = [
        SeqRecord(Seq("ACGT"), id="r1", annotations={}),
        SeqRecord(Seq("ACGT"), id="r2", annotations={"organism": "Influenza A virus (H1N1)"}),
    ]

    assert GenBank._parse_target(records) == "Influenza_A_virus"


def test_split_genbank_creates_fasta_and_gff_with_target(tmp_path: Path) -> None:
    """Test that split_genbank creates FASTA and GFF files with target name.

    Verifies that GenBank.split_genbank() successfully splits a GenBank file
    into separate FASTA and GFF outputs with correct structure and returns
    a normalized target name when emit_target=True.

    Parameters
    ----------
    tmp_path : Path
        Temporary directory provided by pytest.
    """
    source = UNIT_DATA_DIR / "test_reference_double.gb"
    working_copy = tmp_path / source.name
    shutil.copy(source, working_copy)

    fasta_path, gff_path, target = GenBank.split_genbank(working_copy, emit_target=True)

    assert fasta_path == working_copy.with_suffix(".fasta")
    assert gff_path == working_copy.with_suffix(".gff")
    assert target == "Influenza_A_virus"
    assert fasta_path.exists()
    assert gff_path.exists()
    assert fasta_path.stat().st_size > 0
    assert gff_path.stat().st_size > 0

    fasta_records = list(SeqIO.parse(fasta_path, "fasta"))
    assert len(fasta_records) == 2

    gff_text = gff_path.read_text(encoding="utf-8")
    assert "##gff-version" in gff_text


def test_split_genbank_without_target_returns_empty_string(tmp_path: Path) -> None:
    """Test that split_genbank returns empty target when disabled.

    Verifies that GenBank.split_genbank() returns an empty string for the
    target when emit_target=False.

    Parameters
    ----------
    tmp_path : Path
        Temporary directory provided by pytest.
    """
    source = UNIT_DATA_DIR / "test_reference.gb"
    working_copy = tmp_path / source.name
    shutil.copy(source, working_copy)

    _, _, target = GenBank.split_genbank(working_copy, emit_target=False)

    assert target == ""


def test_split_genbank_raises_on_invalid_genbank_file(tmp_path: Path) -> None:
    """Test that split_genbank propagates parser errors.

    Verifies that GenBank.split_genbank() raises ValueError when given a
    malformed GenBank file, propagating parser-derived exceptions with
    descriptive error context.

    Parameters
    ----------
    tmp_path : Path
        Temporary directory provided by pytest.
    """
    source = UNIT_DATA_DIR / "test_reference_bad.gb"
    working_copy = tmp_path / source.name
    shutil.copy(source, working_copy)

    with pytest.raises(ValueError, match=r"Error opening GenBank file"):
        GenBank.split_genbank(working_copy)
