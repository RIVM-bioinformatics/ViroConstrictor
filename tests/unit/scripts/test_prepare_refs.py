import sys
from pathlib import Path

from Bio import SeqIO

from tests.utils.fasta_generator import generate_fasta_file

project_root = Path(__file__).parent.parent.parent.parent
sys.path.insert(0, str(project_root.joinpath("ViroConstrictor/workflow")))
from ViroConstrictor.workflow.main.scripts.prepare_refs import PrepareRefs  # isort:skip


def test_prepare_refs(tmp_path: Path) -> None:
    input_path = tmp_path / "test_ref.fasta"
    output_path = tmp_path / "test_ref_processed.fasta"
    descriptions = ["Test sequence"]
    ids = [x.split(" ", 1)[0] for x in descriptions]  # This is how the reference ID is extracted by SeqIO

    # Create a mock FASTA file
    generate_fasta_file(input_path, length=100, number_of_sequences=1, descriptions=descriptions)

    script = PrepareRefs(input_path, output_path, ids[0])
    script.run()

    output = SeqIO.parse(output_path, "fasta")
    assert output is not None, "Output file should not be empty."
    for record in output:
        assert record.id == ids[0], f"Expected ID {ids[0]}, got {record.id}."
        assert record.seq.isupper(), "Sequence should be in uppercase."
