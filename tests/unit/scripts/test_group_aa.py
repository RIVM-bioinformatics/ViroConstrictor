from pathlib import Path

from ViroConstrictor.workflow.scripts.group_aminoacids import GroupAminoAcids

# def test_group_amino_acids(tmp_path: Path) -> None:
#     input = "tests/unit/data/aa.faa"
#     output = "tests/unit/data/aa_output.faa"
#     pkl = "tests/unit/data/sampleinfo.pkl"

#     a = GroupAminoAcids(input=input, output=output, space=pkl)
#     a.run()


def test_group_amino_acids2(tmp_path: Path) -> None:
    input = "tests/unit/data/aa2.faa"
    result_string = (
        "results/Virus~Severe_acute_respiratory_syndrome_coronavirus_2/RefID~NC_045512.2/aminoacids/ORF7b.faa "
        "results/Virus~Severe_acute_respiratory_syndrome_coronavirus_2/RefID~NC_045512.2/aminoacids/S.faa "
        "results/Virus~Severe_acute_respiratory_syndrome_coronavirus_2/RefID~NC_045512.2/aminoacids/ORF10.faa "
        "results/Virus~Severe_acute_respiratory_syndrome_coronavirus_2/RefID~NC_045512.2/aminoacids/N.faa "
        "results/Virus~Severe_acute_respiratory_syndrome_coronavirus_2/RefID~NC_045512.2/aminoacids/ORF6.faa "
        "results/Virus~Severe_acute_respiratory_syndrome_coronavirus_2/RefID~NC_045512.2/aminoacids/ORF1ab.faa "
        "results/Virus~Severe_acute_respiratory_syndrome_coronavirus_2/RefID~NC_045512.2/aminoacids/ORF3a.faa "
        "results/Virus~Severe_acute_respiratory_syndrome_coronavirus_2/RefID~NC_045512.2/aminoacids/M.faa "
        "results/Virus~Severe_acute_respiratory_syndrome_coronavirus_2/RefID~NC_045512.2/aminoacids/E.faa "
        "results/Virus~Severe_acute_respiratory_syndrome_coronavirus_2/RefID~NC_045512.2/aminoacids/ORF8.faa "
        "results/Virus~Severe_acute_respiratory_syndrome_coronavirus_2/RefID~NC_045512.2/aminoacids/ORF7a.faa"
    )
    pkl = "tests/unit/data/sampleinfo2.pkl"

    a = GroupAminoAcids(input=input, output=result_string, space=pkl)
    a.run()
