from pathlib import Path

from ViroConstrictor.workflow.scripts.group_aminoacids import GroupAminoAcids


def test_group_amino_acids(tmp_path: Path) -> None:
    input = "tests/unit/data/aa.faa"
    output = "tests/unit/data/aa_output.faa"
    pkl = "tests/unit/data/sampleinfo.pkl"

    a = GroupAminoAcids(input=input, output=output, space=pkl)
    a.run()
