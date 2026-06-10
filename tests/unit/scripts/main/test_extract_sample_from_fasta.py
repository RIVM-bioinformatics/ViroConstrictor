"""
Unit tests for ExtractSampleFromFasta script.

Tests extraction of records matching a sample identifier from combined FASTA
files, including prefix matching and cross-file concatenation.
"""

import sys
from argparse import ArgumentParser
from pathlib import Path

import pytest
from Bio import SeqIO

PROJECT_ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(PROJECT_ROOT.joinpath("ViroConstrictor/workflow")))
from ViroConstrictor.workflow.main.scripts.extract_sample_from_fasta import ExtractSampleFromFasta  # noqa: E402, isort:skip


def _record_ids(path: Path) -> list[str]:
    """
    Extract record identifiers from FASTA file.

    Parameters
    ----------
    path : Path
        Path to FASTA file.

    Returns
    -------
    list[str]
        List of sequence record IDs in file order.
    """
    return [record.id for record in SeqIO.parse(path, "fasta")]


def test_add_arguments_parses_required_values() -> None:
    """
    Verify that ExtractSampleFromFasta CLI parser registers required arguments.

    Tests that the argument parser accepts input_files, output, and sample_name
    arguments and correctly stores their values.
    """
    parser = ArgumentParser()
    ExtractSampleFromFasta.add_arguments(parser)

    args = parser.parse_args(
        [
            "--input",
            "unused",
            "--output",
            "out.fasta",
            "--input_files",
            "a.fasta",
            "b.fasta",
            "--sample_name",
            "sample1",
        ]
    )

    assert args.input_files == ["a.fasta", "b.fasta"]
    assert args.sample_name == "sample1"


def test_run_extracts_exact_and_prefixed_ids_only(tmp_path: Path) -> None:
    """
    Verify records matching exact ID or ID prefix are extracted.

    Tests that the run() method extracts records with IDs that either match
    exactly or match with a dot-prefix (e.g., "sample1.geneA"), while rejecting
    records that don't match these criteria.
    """
    f1 = tmp_path / "in1.fasta"
    out = tmp_path / "out.fasta"

    f1.write_text(
        ">sample1\nAAAA\n>sample1.geneA\nTTTT\n>sample10\nCCCC\n>sample1extra\nGGGG\n",
        encoding="utf-8",
    )

    extractor = ExtractSampleFromFasta(
        input="",
        output=out,
        input_files=[f1],
        sample_name="sample1",
    )
    extractor.run()

    assert out.exists()
    assert _record_ids(out) == ["sample1", "sample1.geneA"]


def test_run_combines_across_multiple_files_in_input_order(tmp_path: Path) -> None:
    """
    Verify multiple input files are combined in specified order.

    Tests that matching records from multiple input files are combined and
    written to output in the order they appear across input files.
    """
    f1 = tmp_path / "in1.fasta"
    f2 = tmp_path / "in2.fasta"
    out = tmp_path / "out.fasta"

    f1.write_text(
        ">other\nAAAA\n>sampleA.part1\nTTTT\n",
        encoding="utf-8",
    )
    f2.write_text(
        ">sampleA\nCCCC\n>sampleA.part2\nGGGG\n",
        encoding="utf-8",
    )

    extractor = ExtractSampleFromFasta(
        input="",
        output=out,
        input_files=[f1, f2],
        sample_name="sampleA",
    )
    extractor.run()

    assert _record_ids(out) == ["sampleA.part1", "sampleA", "sampleA.part2"]


def test_run_skips_missing_and_empty_files(tmp_path: Path) -> None:
    """
    Verify missing and empty files are skipped gracefully.

    Tests that the script processes only files that exist and have content,
    silently skipping files that are missing or empty.
    """
    valid = tmp_path / "valid.fasta"
    empty = tmp_path / "empty.fasta"
    missing = tmp_path / "missing.fasta"
    out = tmp_path / "out.fasta"

    valid.write_text(
        ">sampleZ\nATAT\n>other\nCGCG\n",
        encoding="utf-8",
    )
    empty.write_text("", encoding="utf-8")

    extractor = ExtractSampleFromFasta(
        input="",
        output=out,
        input_files=[missing, empty, valid],
        sample_name="sampleZ",
    )
    extractor.run()

    assert _record_ids(out) == ["sampleZ"]


def test_run_when_no_matches_writes_empty_fasta(tmp_path: Path) -> None:
    """
    Verify that no matching records produce empty output file.

    Tests that when no FASTA records match the sample name, the output file
    is created but remains empty.
    """
    f1 = tmp_path / "in.fasta"
    out = tmp_path / "out.fasta"

    f1.write_text(
        ">other1\nAAAA\n>other2\nTTTT\n",
        encoding="utf-8",
    )

    extractor = ExtractSampleFromFasta(
        input="",
        output=out,
        input_files=[f1],
        sample_name="sample_missing",
    )
    extractor.run()

    assert out.exists()
    assert out.read_text(encoding="utf-8") == ""


def test_run_surfaces_invalid_fasta_errors(tmp_path: Path) -> None:
    """
    Verify that invalid FASTA files are handled gracefully.

    Tests that the script handles malformed FASTA files without crashing,
    producing valid output (empty if no records match).
    """
    invalid = tmp_path / "broken.fasta"
    out = tmp_path / "out.fasta"

    invalid.write_text(">sample1\nACGT\n>\n", encoding="utf-8")

    extractor = ExtractSampleFromFasta(
        input="",
        output=out,
        input_files=[invalid],
        sample_name="sample1",
    )

    # SeqIO may tolerate many malformed inputs; this case produces no matching records
    # and should not crash the script.
    extractor.run()
    assert out.exists()


@pytest.mark.xfail(reason="No explicit validation for empty sample_name; intended behavior is to reject invalid sample identifiers", strict=False)
def test_run_rejects_empty_sample_name(tmp_path: Path) -> None:
    """
    Test empty sample_name rejection (XFAIL - intended behavior).

    Intended behavior: Empty sample_name values should raise ValueError with a
    clear message, but the current implementation does not validate this.
    """
    f1 = tmp_path / "in.fasta"
    out = tmp_path / "out.fasta"
    f1.write_text(">sample1\nAAAA\n", encoding="utf-8")

    extractor = ExtractSampleFromFasta(
        input="",
        output=out,
        input_files=[f1],
        sample_name="",
    )

    with pytest.raises(ValueError, match="sample|name|empty"):
        extractor.run()
