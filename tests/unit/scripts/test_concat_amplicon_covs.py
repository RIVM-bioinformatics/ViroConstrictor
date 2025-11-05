import sys
from pathlib import Path

project_root = Path(__file__).parent.parent.parent.parent
sys.path.insert(0, str(project_root.joinpath("ViroConstrictor/workflow")))
from ViroConstrictor.workflow.main.scripts.concat_amplicon_covs import ConcatAmpliconCovs  # isort:skip


def test_amplicon_covs(tmp_path: Path) -> None:

    input = "tests/unit/data/ESIB_EQA_2024_SARS1_01_ampliconcoverage.csv"
    output = "tests/unit/data/ESIB_EQA_2024_SARS1_01_ampliconcoverage_output.csv"

    a = ConcatAmpliconCovs(input=input, output=output)
    a.run()
