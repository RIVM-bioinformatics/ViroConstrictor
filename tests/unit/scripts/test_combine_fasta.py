import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PROJECT_ROOT.joinpath("ViroConstrictor/workflow")))
from ViroConstrictor.workflow.main.scripts.combine_fasta import CombineFasta  # isort:skip

# TODO: These tests need to be verified in more detail as there were made quickly in order to cover the new script.

# TODO: We need to add tests for error handling and invalid inputs (or more broadly unhappy flows).


def test_combine_fasta_basic(tmp_path: Path) -> None:
    """Test basic FASTA combination functionality."""
    input1 = tmp_path / "sample1.fasta"
    input2 = tmp_path / "sample2.fasta"
    output = tmp_path / "combined.fasta"

    input1.write_text(">sample1 mincov=10\nACGT\n")
    input2.write_text(">sample2 mincov=5\nTGCA\n")

    combiner = CombineFasta(
        input="",
        output=output,
        input_files=[input1, input2],
        virus_list=["VirusA", "VirusB"],
        refid_list=["RefA", "RefB"],
    )
    combiner.run()

    assert output.exists()
    content = output.read_text()
    assert "VirusA" in content
    assert "RefA" in content


def test_combine_fasta_empty_input(tmp_path: Path) -> None:
    """Test FASTA combination with empty input files."""
    input1 = tmp_path / "empty.fasta"
    output = tmp_path / "combined.fasta"

    input1.write_text("")

    combiner = CombineFasta(
        input="",
        output=output,
        input_files=[input1],
        virus_list=["VirusA"],
        refid_list=["RefA"],
    )
    combiner.run()

    assert output.exists()


def test_combine_fasta_multiple_sequences(tmp_path: Path) -> None:
    """Test combining FASTA files with multiple sequences."""
    input1 = tmp_path / "multi1.fasta"
    input2 = tmp_path / "multi2.fasta"
    output = tmp_path / "combined.fasta"

    input1.write_text(">seq1 mincov=10\nACGT\n>seq2 mincov=15\nTGCA\n")
    input2.write_text(">seq3 mincov=20\nGGCC\n")

    combiner = CombineFasta(
        input="",
        output=output,
        input_files=[input1, input2],
        virus_list=["VirusA", "VirusB"],
        refid_list=["RefA", "RefB"],
    )
    combiner.run()

    assert output.exists()
    content = output.read_text()
    assert content.count(">") >= 3  # At least 3 sequences
