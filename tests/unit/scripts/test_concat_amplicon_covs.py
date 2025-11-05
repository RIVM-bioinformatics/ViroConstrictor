import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).parent.parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT.joinpath("ViroConstrictor/workflow")))
from ViroConstrictor.workflow.main.scripts.concat_amplicon_covs import ConcatAmpliconCovs  # isort:skip


def test_amplicon_covs(tmp_path: Path) -> None:

    input_coverage = PROJECT_ROOT / "tests" / "unit" / "data" / "ESIB_EQA_2024_SARS1_01_ampliconcoverage.csv"
    output = PROJECT_ROOT / "tests" / "unit" / "data" / "ESIB_EQA_2024_SARS1_01_ampliconcoverage_output.csv"

    a = ConcatAmpliconCovs(input=[""], input_coverages=[input_coverage], output=output)
    a.run()
