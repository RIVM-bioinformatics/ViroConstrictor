"""Unit tests for count_mapped_reads workflow script.

This module tests the CountMappedReads script, which processes BAM alignment files
to extract and aggregate read statistics (mapped counts, mismatches, sequence identity)
for each reference sequence. Tests cover happy paths, edge cases (unmapped reads,
zero references), error handling (missing NM tag), and known limitations.
"""

import sys
import types
from argparse import ArgumentParser
from pathlib import Path
from typing import Dict, Iterable, List

import pandas as pd
import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(PROJECT_ROOT.joinpath("ViroConstrictor/workflow")))
try:
    import pysam  # noqa: F401
except ModuleNotFoundError:
    pysam_stub = types.ModuleType("pysam")
    pysam_stub.AlignmentFile = None
    sys.modules["pysam"] = pysam_stub

from ViroConstrictor.workflow.match_ref.scripts.count_mapped_reads import CountMappedReads  # noqa: E402, isort:skip


def test_count_mapped_reads_add_arguments_accepts_base_inputs() -> None:
    """Verify that CountMappedReads.add_arguments parses base input/output flags."""
    parser = ArgumentParser()
    CountMappedReads.add_arguments(parser)

    args = parser.parse_args(["--input", "in.bam", "--output", "out.csv"])
    assert args.input == "in.bam"
    assert args.output == "out.csv"


class _FakeRead:
    """Mock BAM read for testing CountMappedReads.

    Simulates pysam.AlignmentFile read objects with configurable unmapped status,
    mismatch count (NM tag), and alignment length.
    """

    def __init__(self, *, is_unmapped: bool, nm: int = 0, query_alignment_length: int = 100) -> None:
        """Initialize a mock read.

        Parameters
        ----------
        is_unmapped : bool
            If True, read is marked as unmapped (UNMAPPED flag set).
        nm : int, optional
            Number of mismatches (NM tag value). Default is 0.
        query_alignment_length : int, optional
            Aligned query sequence length. Default is 100.
        """
        self.is_unmapped = is_unmapped
        self._nm = nm
        self.query_alignment_length = query_alignment_length

    def get_tag(self, tag: str) -> int:
        """Retrieve SAM tag value (currently only supports NM tag).

        Parameters
        ----------
        tag : str
            Tag name (e.g., "NM").

        Returns
        -------
        int
            Tag value for "NM" tag.

        Raises
        ------
        KeyError
            If tag is not "NM".
        """
        if tag != "NM":
            raise KeyError(tag)
        return self._nm


class _FakeAlignmentFile:
    """Mock BAM alignment file for testing CountMappedReads.

    Simulates pysam.AlignmentFile with reference-keyed read storage.
    """

    def __init__(self, reads_by_reference: Dict[str, List[_FakeRead]]) -> None:
        """Initialize a mock alignment file.

        Parameters
        ----------
        reads_by_reference : Dict[str, List[_FakeRead]]
            Mapping from reference sequence name to list of reads aligned to it.
        """
        self.references = tuple(reads_by_reference.keys())
        self._reads_by_reference = reads_by_reference

    def fetch(self, ref: str) -> Iterable[_FakeRead]:
        """Fetch reads aligned to a specific reference.

        Parameters
        ----------
        ref : str
            Reference sequence name.

        Returns
        -------
        Iterable[_FakeRead]
            Iterator over reads aligned to the reference.
        """
        return self._reads_by_reference[ref]


def test_count_mapped_reads_happy_path(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """Verify CountMappedReads correctly aggregates statistics for multiple references.

    Tests that mapped read counts, mismatch averages, and sequence identities are
    correctly calculated and written to CSV output.
    """
    reads_by_reference = {
        "refA": [
            _FakeRead(is_unmapped=False, nm=2, query_alignment_length=100),
            _FakeRead(is_unmapped=False, nm=0, query_alignment_length=100),
            _FakeRead(is_unmapped=True, nm=10, query_alignment_length=100),
        ],
        "refB": [
            _FakeRead(is_unmapped=False, nm=5, query_alignment_length=50),
        ],
    }

    monkeypatch.setattr(
        "ViroConstrictor.workflow.match_ref.scripts.count_mapped_reads.pysam.AlignmentFile",
        lambda _path, _mode: _FakeAlignmentFile(reads_by_reference),
    )

    output = tmp_path / "mapped.csv"
    counter = CountMappedReads(input="fake.bam", output=output)

    counter.run()

    df = pd.read_csv(output)
    assert set(df["Reference"]) == {"refA", "refB"}

    row_a = df[df["Reference"] == "refA"].iloc[0]
    assert row_a["Mapped Reads"] == 2
    assert row_a["Avg. Mismatches per Read"] == pytest.approx(1.0)
    assert row_a["Avg. Sequence Identity"] == pytest.approx(0.99)

    row_b = df[df["Reference"] == "refB"].iloc[0]
    assert row_b["Mapped Reads"] == 1
    assert row_b["Avg. Mismatches per Read"] == pytest.approx(5.0)
    assert row_b["Avg. Sequence Identity"] == pytest.approx(0.9)


def test_count_mapped_reads_includes_references_without_mapped_reads(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Verify that references with zero or only unmapped reads are included in output.

    Tests edge case where all reads are unmapped or reference has no reads.
    """
    reads_by_reference = {
        "ref_with_none": [_FakeRead(is_unmapped=True, nm=3, query_alignment_length=100)],
        "ref_with_zero_records": [],
    }

    monkeypatch.setattr(
        "ViroConstrictor.workflow.match_ref.scripts.count_mapped_reads.pysam.AlignmentFile",
        lambda _path, _mode: _FakeAlignmentFile(reads_by_reference),
    )

    output = tmp_path / "mapped_none.csv"
    CountMappedReads(input="fake.bam", output=output).run()

    df = pd.read_csv(output)
    assert set(df["Reference"]) == {"ref_with_none", "ref_with_zero_records"}
    assert (df["Mapped Reads"] == 0).all()
    assert (df["Avg. Mismatches per Read"] == 0).all()
    assert (df["Avg. Sequence Identity"] == 0).all()


def test_count_mapped_reads_missing_nm_tag_raises(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """Verify that a missing NM (mismatch count) tag raises KeyError."""
    class _ReadWithoutNm(_FakeRead):
        def get_tag(self, tag: str) -> int:
            raise KeyError("NM")

    reads_by_reference = {"ref": [_ReadWithoutNm(is_unmapped=False, query_alignment_length=100)]}
    monkeypatch.setattr(
        "ViroConstrictor.workflow.match_ref.scripts.count_mapped_reads.pysam.AlignmentFile",
        lambda _path, _mode: _FakeAlignmentFile(reads_by_reference),
    )

    with pytest.raises(KeyError, match="NM"):
        CountMappedReads(input="fake.bam", output=tmp_path / "out.csv").run()


@pytest.mark.xfail(
    raises=ZeroDivisionError,
    strict=True,
    reason="Intended behavior: zero-length alignments should be handled gracefully (e.g., skipped or treated as 0 identity).",
)
def test_count_mapped_reads_handles_zero_alignment_length(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """Verify that zero-length alignments raise ZeroDivisionError (known limitation)."""
    reads_by_reference = {
        "ref": [
            _FakeRead(is_unmapped=False, nm=0, query_alignment_length=0),
        ]
    }
    monkeypatch.setattr(
        "ViroConstrictor.workflow.match_ref.scripts.count_mapped_reads.pysam.AlignmentFile",
        lambda _path, _mode: _FakeAlignmentFile(reads_by_reference),
    )

    CountMappedReads(input="fake.bam", output=tmp_path / "out.csv").run()
