import sys
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PROJECT_ROOT.joinpath("ViroConstrictor/workflow")))
from ViroConstrictor.workflow.main.scripts.aggregate_combined_files import AggregateCombinedFiles  # isort:skip

# TODO: These tests need to be verified in more detail as there were made quickly in order to cover the new script.

# TODO: We need to add tests for error handling and invalid inputs (or more broadly unhappy flows).


def test_aggregate_mutations(tmp_path: Path) -> None:
    """Test aggregating mutation files."""
    input1 = tmp_path / "mutations1.tsv"
    input2 = tmp_path / "mutations2.tsv"
    output = tmp_path / "all_mutations.tsv"

    input1.write_text("Sample\tVirus\tReference_ID\tPosition\tReference_Base\tVariant_Base\tDepth\nsample1\tVirusA\tRefA\t100\tA\tT\t50\n")
    input2.write_text("Sample\tVirus\tReference_ID\tPosition\tReference_Base\tVariant_Base\tDepth\nsample2\tVirusB\tRefB\t200\tC\tG\t60\n")

    aggregator = AggregateCombinedFiles(
        input="",
        output=output,
        input_files=[input1, input2],
        file_type="mutations",
    )
    aggregator.run()

    assert output.exists()
    df = pd.read_csv(output, sep="\t")
    assert len(df) == 2
    assert "Virus" in df.columns
    assert "Reference_ID" in df.columns


def test_aggregate_coverage(tmp_path: Path) -> None:
    """Test aggregating coverage files."""
    input1 = tmp_path / "coverage1.tsv"
    input2 = tmp_path / "coverage2.tsv"
    output = tmp_path / "all_coverage.tsv"

    input1.write_text("Sample_name\tVirus\tReference_ID\t100\t200\t300\nsample1\tVirusA\tRefA\t10\t20\t30\n")
    input2.write_text("Sample_name\tVirus\tReference_ID\t100\t200\t300\nsample2\tVirusB\tRefB\t15\t25\t35\n")

    aggregator = AggregateCombinedFiles(
        input="",
        output=output,
        input_files=[input1, input2],
        file_type="coverage",
    )
    aggregator.run()

    assert output.exists()
    df = pd.read_csv(output, sep="\t")
    assert len(df) == 2


def test_aggregate_amplicon_coverage(tmp_path: Path) -> None:
    """Test aggregating amplicon coverage files."""
    input1 = tmp_path / "amplicon1.csv"
    input2 = tmp_path / "amplicon2.csv"
    output = tmp_path / "all_amplicon.csv"

    input1.write_text("Virus,Reference_ID,amplicon,mean_cov,median_cov\nVirusA,RefA,amp1,100,95\n")
    input2.write_text("Virus,Reference_ID,amplicon,mean_cov,median_cov\nVirusB,RefB,amp2,150,145\n")

    aggregator = AggregateCombinedFiles(
        input="",
        output=output,
        input_files=[input1, input2],
        file_type="amplicon_coverage",
        separator=",",
    )
    aggregator.run()

    assert output.exists()
    df = pd.read_csv(output, sep=",")
    assert len(df) == 2


def test_aggregate_empty_files(tmp_path: Path) -> None:
    """Test aggregating when all input files are empty."""
    input1 = tmp_path / "empty1.tsv"
    input2 = tmp_path / "empty2.tsv"
    output = tmp_path / "aggregated.tsv"

    input1.write_text("")
    input2.write_text("")

    aggregator = AggregateCombinedFiles(
        input="",
        output=output,
        input_files=[input1, input2],
        file_type="mutations",
    )
    aggregator.run()

    assert output.exists()


def test_aggregate_mixed_empty_and_data(tmp_path: Path) -> None:
    """Test aggregating when some files are empty and some have data."""
    input1 = tmp_path / "empty.tsv"
    input2 = tmp_path / "data.tsv"
    output = tmp_path / "aggregated.tsv"

    input1.write_text("")
    input2.write_text("Sample\tVirus\tReference_ID\tPosition\tReference_Base\tVariant_Base\tDepth\nsample1\tVirusA\tRefA\t100\tA\tT\t50\n")

    aggregator = AggregateCombinedFiles(
        input="",
        output=output,
        input_files=[input1, input2],
        file_type="mutations",
    )
    aggregator.run()

    assert output.exists()
    df = pd.read_csv(output, sep="\t")
    assert len(df) == 1
