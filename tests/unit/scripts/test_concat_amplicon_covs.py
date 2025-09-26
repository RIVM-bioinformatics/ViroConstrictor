from pathlib import Path

from ViroConstrictor.workflow.main.scripts.concat_amplicon_covs import (
    ConcatAmpliconCovs,
)


def test_amplicon_covs(tmp_path: Path) -> None:
    input = "tests/unit/data/ESIB_EQA_2024_SARS1_01_ampliconcoverage.csv"
    output = "tests/unit/data/ESIB_EQA_2024_SARS1_01_ampliconcoverage_output.csv"

    a = ConcatAmpliconCovs(input=input, output=output)
    a.run()
