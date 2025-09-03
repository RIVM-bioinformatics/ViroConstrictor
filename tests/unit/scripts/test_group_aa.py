from pathlib import Path

from ViroConstrictor.workflow.main.scripts.group_aminoacids import GroupAminoAcids

# def test_group_amino_acids(tmp_path: Path) -> None:
#     input = "tests/unit/data/aa.faa"
#     output = "tests/unit/data/aa_output.faa"
#     pkl = "tests/unit/data/sampleinfo.pkl"

#     a = GroupAminoAcids(input=input, output=output, space=pkl)
#     a.run()


def test_group_amino_acids2(tmp_path: Path) -> None:
    input = "tests/unit/data/aa2.faa"
    result_string = (
        "test/unit/data/ORF7b.faa "
        "test/unit/data/S.faa "
        "test/unit/data/ORF10.faa "
        "test/unit/data/N.faa "
        "test/unit/data/ORF6.faa "
        "test/unit/data/ORF1ab.faa "
        "test/unit/data/ORF3a.faa "
        "test/unit/data/M.faa "
        "test/unit/data/E.faa "
        "test/unit/data/ORF8.faa "
        "test/unit/data/ORF7a.faa"
    )
    pkl = "tests/unit/data/sampleinfo2.pkl"

    a = GroupAminoAcids(input=input, output=result_string, space=pkl)
    a.run()
