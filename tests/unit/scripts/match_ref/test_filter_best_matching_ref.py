"""Unit tests for filter_best_matching_ref workflow script.

This module tests the FilterBestMatchingRef script, which selects the best-matching
reference from a set of candidates based on mapped read count and mismatch metrics.
Tests cover happy paths, tie-breaking, empty/invalid inputs, and missing data errors.
"""

from argparse import ArgumentParser
import sys
from pathlib import Path

import pandas as pd
import pytest
from Bio import SeqIO

PROJECT_ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(PROJECT_ROOT.joinpath("ViroConstrictor/workflow")))
from ViroConstrictor.workflow.match_ref.scripts.filter_best_matching_ref import FilterBestMatchingRef  # isort:skip


def test_filter_best_matching_ref_add_arguments_defines_expected_flags() -> None:
    """Verify FilterBestMatchingRef.add_arguments parses input/output and ref flags."""
    parser = ArgumentParser()
    FilterBestMatchingRef.add_arguments(parser)

    args = parser.parse_args(
        ["--input", "counts.csv", "--output", "filtered.csv", "--inputref", "refs.fasta", "--filtref", "best.fasta"]
    )
    assert args.input == "counts.csv"
    assert args.output == "filtered.csv"
    assert args.inputref == "refs.fasta"
    assert args.filtref == "best.fasta"


def _write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    """Write test FASTA records to a file.

    Parameters
    ----------
    path : Path
        Output file path.
    records : list[tuple[str, str]]
        List of (identifier, sequence) tuples.
    """
    with path.open("w", encoding="utf-8") as handle:
        for identifier, sequence in records:
            handle.write(f">{identifier}\n{sequence}\n")


def test_filter_best_matching_ref_happy_path(tmp_path: Path) -> None:
    """Verify FilterBestMatchingRef selects the reference with highest mapped reads.

    Tests that the best reference (highest mapped reads) is selected and written
    to both counts and FASTA output files.
    """
    counts = tmp_path / "counts.csv"
    refs = tmp_path / "refs.fasta"
    out_ref = tmp_path / "best.fasta"
    out_count = tmp_path / "best.csv"

    pd.DataFrame(
        [
            {"Reference": "ref_low", "Mapped Reads": 10, "Avg. Mismatches per Read": 1.0},
            {"Reference": "ref_high", "Mapped Reads": 50, "Avg. Mismatches per Read": 0.5},
        ]
    ).to_csv(counts, index=False)
    _write_fasta(refs, [("ref_low", "AAAA"), ("ref_high", "CCCC")])

    FilterBestMatchingRef(input=counts, inputref=refs, filtref=out_ref, output=out_count).run()

    out_df = pd.read_csv(out_count)
    assert len(out_df) == 1
    assert out_df.iloc[0]["Reference"] == "ref_high"
    assert out_df.iloc[0]["Mapped Reads"] == 50

    records = list(SeqIO.parse(out_ref, "fasta"))
    assert len(records) == 1
    assert records[0].id == "ref_high"
    assert str(records[0].seq) == "CCCC"


def test_filter_best_matching_ref_tie_keeps_first_sorted_entry(tmp_path: Path) -> None:
    """Verify FilterBestMatchingRef breaks ties by selecting first alphabetically."""
    counts = tmp_path / "counts_tie.csv"
    refs = tmp_path / "refs_tie.fasta"
    out_ref = tmp_path / "best_tie.fasta"
    out_count = tmp_path / "best_tie.csv"

    pd.DataFrame(
        [
            {"Reference": "ref_first", "Mapped Reads": 99},
            {"Reference": "ref_second", "Mapped Reads": 99},
        ]
    ).to_csv(counts, index=False)
    _write_fasta(refs, [("ref_first", "AAAA"), ("ref_second", "TTTT")])

    FilterBestMatchingRef(input=counts, inputref=refs, filtref=out_ref, output=out_count).run()

    out_df = pd.read_csv(out_count)
    assert out_df.iloc[0]["Reference"] == "ref_first"
    records = list(SeqIO.parse(out_ref, "fasta"))
    assert len(records) == 1
    assert records[0].id == "ref_first"


def test_filter_best_matching_ref_missing_required_columns_raises(tmp_path: Path) -> None:
    """Verify that missing required columns in counts CSV raises KeyError."""
    counts = tmp_path / "counts_invalid.csv"
    refs = tmp_path / "refs.fasta"

    pd.DataFrame([{"wrong": "ref1", "also_wrong": 12}]).to_csv(counts, index=False)
    _write_fasta(refs, [("ref1", "AAAA")])

    with pytest.raises(KeyError):
        FilterBestMatchingRef(
            input=counts,
            inputref=refs,
            filtref=tmp_path / "best.fasta",
            output=tmp_path / "best.csv",
        ).run()


def test_filter_best_matching_ref_empty_counts_raises(tmp_path: Path) -> None:
    """Verify that empty counts CSV (no data rows) raises IndexError."""
    counts = tmp_path / "counts_empty.csv"
    refs = tmp_path / "refs.fasta"

    pd.DataFrame(columns=["Reference", "Mapped Reads"]).to_csv(counts, index=False)
    _write_fasta(refs, [("ref1", "AAAA")])

    with pytest.raises(IndexError):
        FilterBestMatchingRef(
            input=counts,
            inputref=refs,
            filtref=tmp_path / "best.fasta",
            output=tmp_path / "best.csv",
        ).run()


@pytest.mark.xfail(
    strict=True,
    reason="Intended behavior: selecting a reference that is absent from FASTA should raise a clear error.",
)
def test_filter_best_matching_ref_missing_selected_reference_in_fasta_raises(tmp_path: Path) -> None:
    """Verify that selecting a reference not in FASTA raises ValueError."""
    counts = tmp_path / "counts_missing_ref.csv"
    refs = tmp_path / "refs_missing_ref.fasta"
    out_ref = tmp_path / "best_missing_ref.fasta"

    pd.DataFrame([{"Reference": "not_present", "Mapped Reads": 10}]).to_csv(counts, index=False)
    _write_fasta(refs, [("other_ref", "AAAA")])

    with pytest.raises(ValueError, match="not_present"):
        FilterBestMatchingRef(
            input=counts,
            inputref=refs,
            filtref=out_ref,
            output=tmp_path / "best_missing_ref.csv",
        ).run()
