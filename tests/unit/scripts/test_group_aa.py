import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).parent.parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT.joinpath("ViroConstrictor/workflow")))
from ViroConstrictor.workflow.main.scripts.group_aminoacids import GroupAminoAcids  # isort:skip

# def test_group_amino_acids(tmp_path: Path) -> None:
#     input = "tests/unit/data/aa.faa"
#     output = "tests/unit/data/aa_output.faa"
#     pkl = "tests/unit/data/sampleinfo.pkl"

#     a = GroupAminoAcids(input=input, output=output, space=pkl)
#     a.run()


def test_group_amino_acids2(tmp_path: Path) -> None:
    input = PROJECT_ROOT / "tests" / "unit" / "data" / "aa2.faa"
    result_string = (
        PROJECT_ROOT / "tests" / "unit" / "data" / "ORF7b.faa",
        PROJECT_ROOT / "tests" / "unit" / "data" / "S.faa",
        PROJECT_ROOT / "tests" / "unit" / "data" / "ORF10.faa",
        PROJECT_ROOT / "tests" / "unit" / "data" / "N.faa",
        PROJECT_ROOT / "tests" / "unit" / "data" / "ORF6.faa",
        PROJECT_ROOT / "tests" / "unit" / "data" / "ORF1ab.faa",
        PROJECT_ROOT / "tests" / "unit" / "data" / "ORF3a.faa",
        PROJECT_ROOT / "tests" / "unit" / "data" / "M.faa",
        PROJECT_ROOT / "tests" / "unit" / "data" / "E.faa",
        PROJECT_ROOT / "tests" / "unit" / "data" / "ORF8.faa",
        PROJECT_ROOT / "tests" / "unit" / "data" / "ORF7a.faa",
    )
    pkl = PROJECT_ROOT / "tests" / "unit" / "data" / "sampleinfo2.pkl"

    a = GroupAminoAcids(input=str(input), output=str(result_string), space=str(pkl))
    a.run()
