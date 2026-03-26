import sys
from argparse import ArgumentParser
from pathlib import Path

from Bio import SeqIO
import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(PROJECT_ROOT.joinpath("ViroConstrictor/workflow")))
from ViroConstrictor.workflow.main.scripts.prepare_refs import PrepareRefs  # isort:skip


def test_add_arguments_parses_reference_id() -> None:
    parser = ArgumentParser()
    PrepareRefs.add_arguments(parser)
    args = parser.parse_args(["--input", "in.fasta", "--output", "out.fasta", "--reference_id", "REF_1"])
    assert args.reference_id == "REF_1"


def test_prepare_refs_selects_target_and_uppercases_sequence(tmp_path: Path) -> None:
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
    input_path = tmp_path / "input.fasta"
    output_path = tmp_path / "output.fasta"
    input_path.write_text(
        ">REF_A\nacgtac\n",
        encoding="utf-8",
    )

    PrepareRefs(input_path, output_path, "DOES_NOT_EXIST").run()

    assert (not output_path.exists()) or output_path.read_text(encoding="utf-8") == ""


def test_prepare_refs_rejects_list_input_type() -> None:
    script = PrepareRefs(input=["in.fasta"], output=Path("out.fasta"), reference_id="REF_A")  # type: ignore[arg-type]
    with pytest.raises(AssertionError, match="Input"):
        script.run()


@pytest.mark.xfail(reason="Reference matching currently uses substring logic; intended behavior is exact reference ID matching", strict=True)
def test_prepare_refs_uses_exact_reference_match(tmp_path: Path) -> None:
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
