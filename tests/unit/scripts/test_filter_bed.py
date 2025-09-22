from pathlib import Path

from ViroConstrictor.workflow.main.scripts.filter_bed_input import FilterBedInput


def test_group_amino_acids(tmp_path: Path) -> None:
    input = "tests/unit/data/test_primers.bed"
    output = "tests/unit/data/test_primers_output.bed"
    ref_id = "DQ415319.1"

    a = FilterBedInput(input=input, output=output, reference_id=ref_id)
    a.run()
