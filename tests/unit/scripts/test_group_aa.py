from pathlib import Path

import pytest

from ViroConstrictor.workflow.scripts.group_aminoacids import GroupAminoAcids


def test_group_amino_acids(tmp_path: Path) -> None:
    input = "tests/unit/data/aa.faa"
    output = "results/Virus~Influenza_A/RefID~DQ415318.1/aminoacids/Taiwan.faa"
    pkl = "tests/unit/data/sampleinfo.pkl"

    a = GroupAminoAcids(input=input, output=output, space=pkl)
    a.run()


pytest.main(["-v", "-s", __file__])
