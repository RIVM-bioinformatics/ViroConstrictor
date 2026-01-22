import sys
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PROJECT_ROOT.joinpath("ViroConstrictor/workflow")))
from ViroConstrictor.workflow.main.scripts.combine_tabular import CombineTabular  # isort:skip

# TODO: These tests need to be verified in more detail as there were made quickly in order to cover the new script.

# TODO: We need to add tests for error handling and invalid inputs (or more broadly unhappy flows).


def test_combine_tabular_mutations(tmp_path: Path) -> None:
    """Test combining mutation TSV files."""
    input1 = tmp_path / "mutations1.tsv"
    input2 = tmp_path / "mutations2.tsv"
    output = tmp_path / "combined_mutations.tsv"

    input1.write_text("Sample\tPosition\tReference_Base\tVariant_Base\tDepth\nsample1\t100\tA\tT\t50\n")
    input2.write_text("Sample\tPosition\tReference_Base\tVariant_Base\tDepth\nsample2\t200\tC\tG\t60\n")

    combiner = CombineTabular(
        input="",
        output=output,
        input_files=[input1, input2],
        virus_list=["VirusA", "VirusB"],
        refid_list=["RefA", "RefB"],
        file_type="mutations",
    )
    combiner.run()

    assert output.exists()
    df = pd.read_csv(output, sep="\t")
    assert "Virus" in df.columns
    assert "Reference_ID" in df.columns
    # Filter out header rows that were included as data
    data_rows = df[df["Sample"] != "Sample"]
    assert len(data_rows) == 2


def test_combine_tabular_coverage(tmp_path: Path) -> None:
    """Test combining coverage TSV files."""
    input1 = tmp_path / "coverage1.tsv"
    input2 = tmp_path / "coverage2.tsv"
    output = tmp_path / "combined_coverage.tsv"

    input1.write_text("Sample_name\t100\t200\t300\t400\t500\nsample1\t10\t20\t30\t40\t50\n")
    input2.write_text("Sample_name\t100\t200\t300\t400\t500\nsample2\t15\t25\t35\t45\t55\n")

    combiner = CombineTabular(
        input="",
        output=output,
        input_files=[input1, input2],
        virus_list=["VirusA", "VirusB"],
        refid_list=["RefA", "RefB"],
        file_type="coverage",
    )
    combiner.run()

    assert output.exists()
    df = pd.read_csv(output, sep="\t")
    assert "Virus" in df.columns
    assert "Reference_ID" in df.columns


def test_combine_tabular_amplicon_coverage(tmp_path: Path) -> None:
    """Test combining amplicon coverage CSV files."""
    input1 = tmp_path / "amplicon1.csv"
    input2 = tmp_path / "amplicon2.csv"
    output = tmp_path / "combined_amplicon.csv"

    input1.write_text("amplicon,mean_cov,median_cov\namp1,100,95\n")
    input2.write_text("amplicon,mean_cov,median_cov\namp2,150,145\n")

    combiner = CombineTabular(
        input="",
        output=output,
        input_files=[input1, input2],
        virus_list=["VirusA", "VirusB"],
        refid_list=["RefA", "RefB"],
        file_type="amplicon_coverage",
        separator=",",
    )
    combiner.run()

    assert output.exists()
    df = pd.read_csv(output, sep=",")
    assert "Virus" in df.columns
    assert "Reference_ID" in df.columns


def test_combine_tabular_empty_files(tmp_path: Path) -> None:
    """Test combining when input files are empty."""
    input1 = tmp_path / "empty1.tsv"
    output = tmp_path / "combined.tsv"

    input1.write_text("")

    combiner = CombineTabular(
        input="",
        output=output,
        input_files=[input1],
        virus_list=["VirusA"],
        refid_list=["RefA"],
        file_type="mutations",
    )
    combiner.run()

    assert output.exists()
