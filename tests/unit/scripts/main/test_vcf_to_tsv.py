import sys
from argparse import ArgumentParser
from pathlib import Path

import pandas as pd
import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(PROJECT_ROOT.joinpath("ViroConstrictor/workflow")))
from ViroConstrictor.workflow.main.scripts.vcf_to_tsv import VcfToTsv  # isort:skip


def test_add_arguments_parses_required_fields() -> None:
    parser = ArgumentParser()
    VcfToTsv.add_arguments(parser)
    args = parser.parse_args(["--input", "in.vcf", "--output", "out.tsv", "--samplename", "S1"])
    assert args.samplename == "S1"


def test_run_converts_variants_and_filters_unknown_alt(tmp_path: Path) -> None:
    input_path = tmp_path / "example.vcf"
    output_path = tmp_path / "result.tsv"

    input_path.write_text(
        """##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t12345\t.\tA\tT\t.\tPASS\tDP=100;AF=0.90
chr2\t67890\t.\tG\tN\t.\tPASS\tDP=50
""",
        encoding="utf-8",
    )

    VcfToTsv(input_path, output_path, "Sample1").run()

    df = pd.read_csv(output_path, sep="\t", header=None)
    assert df.shape == (1, 6)
    assert df.iloc[0].tolist() == ["Sample1", "chr1", 12345, "A", "T", 100]


def test_run_writes_empty_output_when_all_alt_are_filtered(tmp_path: Path) -> None:
    input_path = tmp_path / "only_unknown.vcf"
    output_path = tmp_path / "result.tsv"
    input_path.write_text(
        """##fileformat=VCFv4.2
chr1\t12345\t.\tA\tN\t.\tPASS\tDP=100
""",
        encoding="utf-8",
    )

    VcfToTsv(input_path, output_path, "Sample1").run()
    assert output_path.exists()
    assert output_path.read_text(encoding="utf-8") == ""


def test_run_with_no_variant_rows_creates_empty_output(tmp_path: Path) -> None:
    input_path = tmp_path / "headers_only.vcf"
    output_path = tmp_path / "result.tsv"
    input_path.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n", encoding="utf-8")

    VcfToTsv(input_path, output_path, "Sample1").run()
    assert output_path.exists()
    assert output_path.read_text(encoding="utf-8") == ""


@pytest.mark.xfail(reason="Malformed INFO fields currently raise during split; intended behavior is robust parsing or clear validation error", strict=False)
def test_run_handles_malformed_info_field_with_clear_error(tmp_path: Path) -> None:
    input_path = tmp_path / "bad_info.vcf"
    output_path = tmp_path / "result.tsv"
    input_path.write_text(
        "chr1\t12345\t.\tA\tT\t.\tPASS\tMALFORMED_INFO\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="INFO|format|malformed"):
        VcfToTsv(input_path, output_path, "Sample1").run()
