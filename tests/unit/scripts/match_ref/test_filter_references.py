"""Unit tests for filter_references workflow script.

This module tests the FilterReferences script, which filters FASTA records by
segment metadata embedded in sequence descriptions. Tests cover segment-based
filtering, "None" wildcard handling, validation errors, and malformed records.
"""

import sys
from argparse import ArgumentParser
from pathlib import Path

import pytest
from Bio import SeqIO

PROJECT_ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(PROJECT_ROOT.joinpath("ViroConstrictor/workflow")))
from ViroConstrictor.workflow.match_ref.scripts.filter_references import FilterReferences  # noqa: E402, isort:skip


def _write_fasta(path: Path, records: list[tuple[str, str, str]]) -> None:
    """Write test FASTA records to a file.

    Parameters
    ----------
    path : Path
        Output file path.
    records : list[tuple[str, str, str]]
        List of (identifier, description, sequence) tuples.
    """
    with path.open("w", encoding="utf-8") as handle:
        for identifier, description, sequence in records:
            handle.write(f">{identifier} {description}\n{sequence}\n")


def test_add_arguments_parses_wildcard_segment() -> None:
    """Verify FilterReferences.add_arguments parses input/output and wildcard_segment."""
    parser = ArgumentParser()
    FilterReferences.add_arguments(parser)

    args = parser.parse_args(["--input", "in.fasta", "--output", "out.fasta", "--wildcard_segment", "HA"])
    assert args.input == "in.fasta"
    assert args.output == "out.fasta"
    assert args.wildcard_segment == "HA"


def test_run_filters_records_by_segment(tmp_path: Path) -> None:
    """Verify FilterReferences filters records by segment metadata in description.

    Tests that only sequences with matching segment (HA/NA) in description are kept.
    """
    input_fasta = tmp_path / "references.fasta"
    output_fasta = tmp_path / "filtered.fasta"

    _write_fasta(
        input_fasta,
        [
            ("ref1", "HA|strainA", "AAAA"),
            ("ref2", "NA|strainB", "CCCC"),
            ("ref3", "HA|strainC", "GGGG"),
        ],
    )

    FilterReferences(input=input_fasta, output=output_fasta, wildcard_segment="HA").run()

    kept = list(SeqIO.parse(output_fasta, "fasta"))
    assert [record.id for record in kept] == ["ref1", "ref3"]


def test_run_with_none_wildcard_keeps_all_records(tmp_path: Path) -> None:
    """Verify FilterReferences with 'None' wildcard keeps all records unfiltered."""
    input_fasta = tmp_path / "references.fasta"
    output_fasta = tmp_path / "all.fasta"

    _write_fasta(
        input_fasta,
        [
            ("ref1", "HA|strainA", "AAAA"),
            ("ref2", "NA|strainB", "CCCC"),
        ],
    )

    FilterReferences(input=str(input_fasta), output=str(output_fasta), wildcard_segment="None").run()

    kept = list(SeqIO.parse(output_fasta, "fasta"))
    assert len(kept) == 2


def test_run_with_missing_input_raises_file_not_found(tmp_path: Path) -> None:
    """Verify that a missing input FASTA file raises FileNotFoundError."""
    with pytest.raises(FileNotFoundError):
        FilterReferences(
            input=tmp_path / "does_not_exist.fasta",
            output=tmp_path / "unused.fasta",
            wildcard_segment="HA",
        ).run()


def test_run_rejects_non_string_wildcard_type(tmp_path: Path) -> None:
    """Verify that non-string wildcard_segment value raises AssertionError."""
    input_fasta = tmp_path / "references.fasta"
    output_fasta = tmp_path / "filtered.fasta"
    _write_fasta(input_fasta, [("ref1", "HA|strainA", "AAAA")])

    script = FilterReferences(input=input_fasta, output=output_fasta, wildcard_segment="HA")
    script.wildcard_segment = 42  # type: ignore[assignment]

    with pytest.raises(AssertionError, match="Wildcard segment"):
        script.run()


@pytest.mark.xfail(
    strict=True,
    reason="Intended behavior: records lacking segment metadata should be skipped, not crash filtering.",
)
def test_run_skips_records_with_malformed_description_when_filtering(tmp_path: Path) -> None:
    """Verify that malformed records (missing segment) are skipped during filtering."""
    input_fasta = tmp_path / "references_malformed.fasta"
    output_fasta = tmp_path / "filtered.fasta"

    with input_fasta.open("w", encoding="utf-8") as handle:
        handle.write(">ref1\nAAAA\n")
        handle.write(">ref2 HA|strainB\nCCCC\n")

    FilterReferences(input=input_fasta, output=output_fasta, wildcard_segment="HA").run()

    kept = list(SeqIO.parse(output_fasta, "fasta"))
    assert [record.id for record in kept] == ["ref2"]
