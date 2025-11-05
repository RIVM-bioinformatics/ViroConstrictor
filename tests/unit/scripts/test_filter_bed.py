import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).parent.parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT.joinpath("ViroConstrictor/workflow")))
from ViroConstrictor.workflow.main.scripts.filter_bed_input import FilterBedInput  # isort:skip


def test_group_amino_acids(tmp_path: Path) -> None:
    input = PROJECT_ROOT / "tests" / "unit" / "data" / "test_primers.bed"
    output = PROJECT_ROOT / "tests" / "unit" / "data" / "test_primers_output.bed"
    ref_id = "DQ415319.1"

    a = FilterBedInput(input=str(input), output=str(output), reference_id=ref_id)
    a.run()
