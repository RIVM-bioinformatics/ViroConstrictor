"""
Unit tests for PrepareRefs script.

Tests reference sequence preparation: selection by ID, uppercase conversion,
and handling of missing/non-matching references.
"""

import sys
from argparse import ArgumentParser
from pathlib import Path

import pytest
from Bio import SeqIO

PROJECT_ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(PROJECT_ROOT.joinpath("ViroConstrictor/workflow")))
from ViroConstrictor.workflow.main.scripts.prepare_refs import PrepareRefs  # isort:skip


def test_add_arguments_parses_reference_id() -> None:
    """
    Verify that PrepareRefs CLI parser registers reference_id argument.

    Tests that the argument parser accepts input, output, and reference_id
    arguments and correctly stores their values.
    """
    parser = ArgumentParser()
    PrepareRefs.add_arguments(parser)
    args = parser.parse_args(["--input", "in.fasta", "--output", "out.fasta", "--reference_id", "REF_1"])
    assert args.reference_id == "REF_1"


def test_prepare_refs_selects_target_and_uppercases_sequence(tmp_path: Path) -> None:
    """
    Verify target reference is selected and sequence converted to uppercase.

    Tests that the run() method extracts the specified reference sequence
    and converts all nucleotides to uppercase.
    """
    input_path = tmp_path / "input.fasta"
    output_path = tmp_path / "output.fasta"
    input_path.write_text(
        ">REF_A\nacgtac\n>REF_B\nttttgg\n",
        encoding="utf-8",
    )

    PrepareRefs(input_path, output_path, "REF_A").run()

    records = list(SeqIO.parse(output_path, "fasta"))
    assert len(records) == 1
    assert records[0].id == "REF_A"
    assert str(records[0].seq) == "ACGTAC"


def test_prepare_refs_no_match_does_not_write_sequence(tmp_path: Path) -> None:
    """
    Verify that non-matching reference produces no output.

    Tests that when no reference sequence matches the specified reference_id,
    the output file is not created or remains empty.
    """
    input_path = tmp_path / "input.fasta"
    output_path = tmp_path / "output.fasta"
    input_path.write_text(
        ">REF_A\nacgtac\n",
        encoding="utf-8",
    )

    PrepareRefs(input_path, output_path, "DOES_NOT_EXIST").run()

    assert (not output_path.exists()) or output_path.read_text(encoding="utf-8") == ""


def test_prepare_refs_rejects_list_input_type() -> None:
    """
    Verify that list input type raises AssertionError.

    Tests that the run() method validates input is a string or Path, not a list,
    and raises AssertionError for type mismatches.
    """
    script = PrepareRefs(input=["in.fasta"], output=Path("out.fasta"), reference_id="REF_A")  # type: ignore[arg-type]
    with pytest.raises(AssertionError, match="Input"):
        script.run()


@pytest.mark.xfail(reason="Reference matching currently uses substring logic; intended behavior is exact reference ID matching", strict=True)
def test_prepare_refs_uses_exact_reference_match(tmp_path: Path) -> None:
    """
    Test exact reference matching (XFAIL - intended behavior, strict=True).

    Intended behavior: Reference matching should be exact (REF_1 != REF_10),
    but the current implementation uses substring matching which is too permissive.
    """
    input_path = tmp_path / "input.fasta"
    output_path = tmp_path / "output.fasta"
    input_path.write_text(
        ">REF_1\nacgt\n>REF_10\ntttt\n",
        encoding="utf-8",
    )

    PrepareRefs(input_path, output_path, "REF_1").run()
    records = list(SeqIO.parse(output_path, "fasta"))
    assert len(records) == 1
    assert records[0].id == "REF_1"
