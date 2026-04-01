"""
Unit tests for Clipper script.

Tests BAM read filtering, trimming, region filtering, spliced read detection,
and CIGAR string parsing for quality-based read filtering workflows.
"""

import sys
import types
from argparse import ArgumentParser
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(PROJECT_ROOT.joinpath("ViroConstrictor/workflow")))

# Keep unit tests hermetic when pysam is unavailable in the active environment.
if "pysam" not in sys.modules:
    sys.modules["pysam"] = types.SimpleNamespace(AlignmentFile=None)

from ViroConstrictor.workflow.main.scripts.clipper import Clipper  # noqa: E402, isort:skip


class FakeRead:
    """
    Mock SAM read record for testing.

    Simulates essential attributes of a pysam AlignmentFile read to allow
    Clipper to be tested without requiring actual BAM file I/O.
    """

    def __init__(
        self,
        query_name: str,
        sequence: str,
        quality: str,
        query_alignment_start: int,
        query_alignment_end: int,
        reference_start: int,
        reference_end: int,
        query_alignment_length: int,
        cigarstring: str = "",
        is_reverse: bool = False,
    ) -> None:
        """
        Initialize a FakeRead with read attributes.

        Parameters
        ----------
        query_name : str
            Read sequence identifier.
        sequence : str
            Aligned query sequence.
        quality : str
            Quality scores string.
        query_alignment_start : int
            Start position on query sequence (0-based).
        query_alignment_end : int
            End position on query sequence (0-based, exclusive).
        reference_start : int
            Start position on reference sequence (0-based).
        reference_end : int
            End position on reference sequence (0-based, exclusive).
        query_alignment_length : int
            Length of query alignment.
        cigarstring : str, optional
            CIGAR string (e.g., "10M3I2D"). Default: "".
        is_reverse : bool, optional
            Whether read is reverse-complemented. Default: False.
        """
        self.query_name = query_name
        self.query_alignment_sequence = sequence
        self.qual = quality
        self.query_alignment_start = query_alignment_start
        self.query_alignment_end = query_alignment_end
        self.reference_start = reference_start
        self.reference_end = reference_end
        self.query_alignment_length = query_alignment_length
        self.cigarstring = cigarstring
        self.is_reverse = is_reverse


class FakeAlignmentFile:
    """
    Mock BAM/SAM alignment file for testing.

    Simulates pysam.AlignmentFile interface to allow testing without real BAM files.
    """

    def __init__(self, reads: list[FakeRead], reflength: int = 1000) -> None:
        """
        Initialize with a collection of FakeRead records.

        Parameters
        ----------
        reads : list[FakeRead]
            List of read records to iterate over.
        reflength : int, optional
            Reference sequence length. Default: 1000.
        """
        self._reads = reads
        self.lengths = [reflength]

    def __iter__(self):
        """Iterate over contained reads."""
        return iter(self._reads)


def _parse_fastq(output_path: Path) -> list[tuple[str, str, str]]:
    """
    Parse FASTQ output file into (name, sequence, quality) tuples.

    Parameters
    ----------
    output_path : Path
        Path to FASTQ file.

    Returns
    -------
    list[tuple[str, str, str]]
        List of (read_name, sequence, quality) tuples.
    """
    lines = output_path.read_text(encoding="utf-8").splitlines()
    records = []
    for i in range(0, len(lines), 4):
        name = lines[i][1:]
        seq = lines[i + 1]
        qual = lines[i + 3]
        records.append((name, seq, qual))
    return records


def _mock_alignment_file(monkeypatch: pytest.MonkeyPatch, reads: list[FakeRead], reflength: int = 1000) -> None:
    """
    Monkey-patch pysam.AlignmentFile to return FakeAlignmentFile.

    Parameters
    ----------
    monkeypatch : pytest.MonkeyPatch
        Pytest fixture for patching.
    reads : list[FakeRead]
        Reads to include in mocked file.
    reflength : int, optional
        Reference length to simulate. Default: 1000.
    """

    def _fake_alignment_file(input_path: Path | str, mode: str, threads: int = 1) -> FakeAlignmentFile:
        assert mode == "rb"
        assert threads >= 1
        assert str(input_path)
        return FakeAlignmentFile(reads=reads, reflength=reflength)

    monkeypatch.setattr("ViroConstrictor.workflow.main.scripts.clipper.pysam.AlignmentFile", _fake_alignment_file)


def test_add_arguments_parses_all_supported_flags() -> None:
    """
    Verify that Clipper CLI parser registers all filtering and region flags.

    Tests that the argument parser accepts input, output, read length filtering,
    spliced read filtering, region filtering, and threading flags.
    """
    parser = ArgumentParser()
    Clipper.add_arguments(parser)

    args = parser.parse_args(
        [
            "--input",
            "in.bam",
            "--output",
            "out.fastq",
            "--exclude-spliced",
            "--spliced-length-threshold",
            "30",
            "--min-aligned-length",
            "0.2",
            "--max-aligned-length",
            "500",
            "--only-include-region",
            "100:900",
            "--threads",
            "4",
        ]
    )

    assert args.input == "in.bam"
    assert args.output == "out.fastq"
    assert args.exclude_spliced is True
    assert args.spliced_length_threshold == 30
    assert args.min_aligned_length == 0.2
    assert args.max_aligned_length == 500.0
    assert args.only_include_region == "100:900"
    assert args.threads == 4


def test_split_cigar_parses_multidigit_operations() -> None:
    """
    Verify CIGAR string parsing correctly handles multi-digit operation counts.

    Tests that the parser breaks down a CIGAR string into (count, operation) tuples,
    including operations with counts > 9.
    """
    assert Clipper._split_cigar("10M3I2D45N6S") == [(10, "M"), (3, "I"), (2, "D"), (45, "N"), (6, "S")]


def test_is_spliced_detects_presence_of_n_operation() -> None:
    """
    Verify spliced read detection identifies N (skip) operations in CIGAR.

    Tests that reads containing N (intron skip) operations are correctly identified
    as spliced, while other operations do not trigger spliced detection.
    """
    assert Clipper._is_spliced([(10, "M"), (20, "N"), (5, "M")]) is True
    assert Clipper._is_spliced([(10, "M"), (2, "I"), (5, "D")]) is False


def test_get_largest_spliced_len_returns_largest_internal_block() -> None:
    """
    Verify largest spliced region block size is correctly identified.

    Tests that the function identifies and returns the largest N (skip) operation
    count from a CIGAR tuple list.
    """
    cigar_tuples = [(10, "M"), (20, "N"), (5, "M"), (50, "N"), (10, "M")]
    assert Clipper._get_largest_spliced_len(cigar_tuples) == 50


@pytest.mark.xfail(reason="Trailing N block is not accounted for in largest-splice calculation", strict=False)
def test_get_largest_spliced_len_counts_trailing_splice_block() -> None:
    """
    Test trailing splice block handling (XFAIL - intended behavior).

    Intended behavior: Trailing spliced region (N operations at the end) should
    contribute to the max splice size calculation.
    """
    # Intended behavior: trailing spliced region should contribute to max splice size.
    cigar_tuples = [(10, "M"), (30, "N")]
    assert Clipper._get_largest_spliced_len(cigar_tuples) == 30


def test_run_writes_forward_and_reverse_reads_with_expected_orientation(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """
    Verify forward and reverse reads are written with correct orientation and reverse-complement.

    Tests that forward reads are preserved as-is while reverse reads are
    reverse-complemented in the output FASTQ.
    """
    output = tmp_path / "reads.fastq"
    reads = [
        FakeRead(
            query_name="read-forward",
            sequence="AACCGG",
            quality="abcdef",
            query_alignment_start=0,
            query_alignment_end=6,
            reference_start=100,
            reference_end=106,
            query_alignment_length=6,
            cigarstring="6M",
        ),
        FakeRead(
            query_name="read-reverse",
            sequence="AGTN",
            quality="1234",
            query_alignment_start=0,
            query_alignment_end=4,
            reference_start=120,
            reference_end=124,
            query_alignment_length=4,
            cigarstring="4M",
            is_reverse=True,
        ),
    ]
    _mock_alignment_file(monkeypatch, reads=reads, reflength=100)

    Clipper(input=tmp_path / "in.bam", output=output).run()

    fastq_records = _parse_fastq(output)
    assert fastq_records[0] == ("read-forward", "AACCGG", "abcdef")
    assert fastq_records[1] == ("read-reverse", "NACT", "4321")


def test_run_filters_empty_and_spliced_reads_when_requested(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """
    Verify empty and spliced reads are filtered when exclude_spliced flag is set.

    Tests that reads with zero alignment length and reads with large spliced
    regions are removed from the output when filtering is enabled.
    """
    output = tmp_path / "reads.fastq"
    reads = [
        FakeRead(
            query_name="empty-trimmed",
            sequence="",
            quality="",
            query_alignment_start=0,
            query_alignment_end=0,
            reference_start=50,
            reference_end=50,
            query_alignment_length=0,
            cigarstring="0M",
        ),
        FakeRead(
            query_name="spliced-large",
            sequence="A" * 40,
            quality="I" * 40,
            query_alignment_start=0,
            query_alignment_end=40,
            reference_start=60,
            reference_end=100,
            query_alignment_length=40,
            cigarstring="20M70N20M",
        ),
        FakeRead(
            query_name="kept",
            sequence="C" * 40,
            quality="J" * 40,
            query_alignment_start=0,
            query_alignment_end=40,
            reference_start=70,
            reference_end=110,
            query_alignment_length=40,
            cigarstring="40M",
        ),
    ]
    _mock_alignment_file(monkeypatch, reads=reads, reflength=500)

    Clipper(input=tmp_path / "in.bam", output=output, exclude_spliced=True, spliced_length_threshold=50).run()

    fastq_records = _parse_fastq(output)
    assert [name for name, _, _ in fastq_records] == ["kept"]


def test_run_applies_region_filter(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """
    Verify region filter restricts output to specified reference interval.

    Tests that only reads whose alignment overlaps the specified genomic region
    are included in the output.
    """
    output = tmp_path / "reads.fastq"
    reads = [
        FakeRead(
            query_name="in-region",
            sequence="ACGT",
            quality="!!!!",
            query_alignment_start=0,
            query_alignment_end=4,
            reference_start=100,
            reference_end=104,
            query_alignment_length=4,
            cigarstring="4M",
        ),
        FakeRead(
            query_name="left-outside",
            sequence="TTTT",
            quality="####",
            query_alignment_start=0,
            query_alignment_end=4,
            reference_start=99,
            reference_end=103,
            query_alignment_length=4,
            cigarstring="4M",
        ),
        FakeRead(
            query_name="right-outside",
            sequence="GGGG",
            quality="$$$$",
            query_alignment_start=0,
            query_alignment_end=4,
            reference_start=120,
            reference_end=151,
            query_alignment_length=4,
            cigarstring="4M",
        ),
    ]
    _mock_alignment_file(monkeypatch, reads=reads, reflength=1000)

    Clipper(input=tmp_path / "in.bam", output=output, only_include_region="100:150").run()

    fastq_records = _parse_fastq(output)
    assert [name for name, _, _ in fastq_records] == ["in-region"]


@pytest.mark.xfail(reason="Length threshold checks are exclusive; intended behavior is inclusive", strict=False)
def test_run_includes_reads_equal_to_min_and_max_length_thresholds(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """
    Test inclusive read length threshold filtering (XFAIL - intended behavior).

    Intended behavior: Reads with lengths exactly equal to minimum and maximum
    thresholds should be retained, but the implementation uses exclusive checks.
    """
    # Intended behavior: reads exactly at minimum and maximum length should be retained.
    output = tmp_path / "reads.fastq"
    reads = [
        FakeRead(
            query_name="equal-min",
            sequence="A" * 100,
            quality="I" * 100,
            query_alignment_start=0,
            query_alignment_end=100,
            reference_start=100,
            reference_end=200,
            query_alignment_length=100,
            cigarstring="100M",
        ),
        FakeRead(
            query_name="equal-max",
            sequence="C" * 200,
            quality="J" * 200,
            query_alignment_start=0,
            query_alignment_end=200,
            reference_start=200,
            reference_end=400,
            query_alignment_length=200,
            cigarstring="200M",
        ),
    ]
    _mock_alignment_file(monkeypatch, reads=reads, reflength=1000)

    Clipper(input=tmp_path / "in.bam", output=output, min_aligned_length=0.1, max_aligned_length=200).run()

    fastq_records = _parse_fastq(output)
    assert [name for name, _, _ in fastq_records] == ["equal-min", "equal-max"]


@pytest.mark.xfail(reason="Malformed region format raises IndexError instead of clear ValueError", strict=False)
def test_run_rejects_malformed_region_with_clear_error(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """
    Test malformed region format rejection (XFAIL - intended behavior).

    Intended behavior: Input validation should raise a clear ValueError for
    malformed region format, but currently raises IndexError instead.
    """
    # Intended behavior: input validation should raise a clear ValueError for malformed regions.
    output = tmp_path / "reads.fastq"
    reads = [
        FakeRead(
            query_name="any",
            sequence="ACGT",
            quality="IIII",
            query_alignment_start=0,
            query_alignment_end=4,
            reference_start=100,
            reference_end=104,
            query_alignment_length=4,
            cigarstring="4M",
        )
    ]
    _mock_alignment_file(monkeypatch, reads=reads, reflength=1000)

    with pytest.raises(ValueError, match="region|format|start:end"):
        Clipper(input=tmp_path / "in.bam", output=output, only_include_region="100-200").run()
