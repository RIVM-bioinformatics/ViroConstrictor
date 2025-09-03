from pathlib import Path

from ViroConstrictor.workflow.scripts.amplicon_covs import AmpliconCovs


def test_amplicon_covs(tmp_path: Path) -> None:
    input = "tests/unit/data/ESIB_EQA_2024_SARS1_01_primers.bed"
    output = "tests/unit/data/ESIB_EQA_2024_SARS1_01_primers_output.csv"
    coverages = "tests/unit/data/ESIB_EQA_2024_SARS1_01_coverage.tsv"
    key = "ESIB_EQA_2024_SARS1_01"

    a = AmpliconCovs(input=input, output=output, key=key, coverages=coverages)
    a.run()
