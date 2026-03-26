import sys
from argparse import ArgumentParser
from pathlib import Path

import pytest
from Bio import SeqIO

PROJECT_ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(PROJECT_ROOT.joinpath("ViroConstrictor/workflow")))
from ViroConstrictor.workflow.main.scripts.combine_fasta import CombineFasta  # isort:skip

def _parse_fasta(path: Path) -> list[tuple[str, str, str]]:
    records = []
    for record in SeqIO.parse(path, "fasta"):
        records.append((record.id, record.description, str(record.seq)))
    return records


def test_add_arguments_parses_required_lists() -> None:
    parser = ArgumentParser()
    CombineFasta.add_arguments(parser)

    args = parser.parse_args(
        [
            "--input",
            "unused",
            "--output",
            "combined.fasta",
            "--input_files",
            "a.fasta",
            "b.fasta",
            "--virus_list",
            "VirusA",
            "VirusB",
            "--refid_list",
            "RefA",
            "RefB",
        ]
    )

    assert args.input_files == ["a.fasta", "b.fasta"]
    assert args.virus_list == ["VirusA", "VirusB"]
    assert args.refid_list == ["RefA", "RefB"]


def test_run_combines_and_rewrites_headers_for_mincov_and_non_mincov(tmp_path: Path) -> None:
    input1 = tmp_path / "first.fasta"
    input2 = tmp_path / "second.fasta"
    output = tmp_path / "combined.fasta"

    input1.write_text(">sample1 mincov=10\nACGT\n>sample2 source=labA\nTTAA\n", encoding="utf-8")
    input2.write_text(">sample3\nGGCC\n", encoding="utf-8")

    combiner = CombineFasta(
        input="",
        output=output,
        input_files=[input1, input2],
        virus_list=["VirusA", "VirusB"],
        refid_list=["RefA", "RefB"],
    )
    combiner.run()

    records = _parse_fasta(output)
    assert [rec[0] for rec in records] == ["sample1", "sample2", "sample3"]
    assert "sample1 VirusA RefA mincov=10" in records[0][1]
    assert "sample2 source=labA VirusA RefA" in records[1][1]
    assert "sample3 VirusB RefB" in records[2][1]
    assert [rec[2] for rec in records] == ["ACGT", "TTAA", "GGCC"]


def test_run_skips_missing_and_empty_inputs(tmp_path: Path) -> None:
    valid = tmp_path / "valid.fasta"
    empty = tmp_path / "empty.fasta"
    missing = tmp_path / "missing.fasta"
    output = tmp_path / "combined.fasta"

    valid.write_text(">sampleA mincov=7\nACAC\n", encoding="utf-8")
    empty.write_text("", encoding="utf-8")

    combiner = CombineFasta(
        input="",
        output=output,
        input_files=[missing, empty, valid],
        virus_list=["VirusM", "VirusE", "VirusV"],
        refid_list=["RefM", "RefE", "RefV"],
    )
    combiner.run()

    records = _parse_fasta(output)
    assert len(records) == 1
    assert records[0][0] == "sampleA"
    assert "VirusV RefV mincov=7" in records[0][1]


def test_run_with_only_missing_inputs_produces_empty_output(tmp_path: Path) -> None:
    missing1 = tmp_path / "missing1.fasta"
    missing2 = tmp_path / "missing2.fasta"
    output = tmp_path / "combined.fasta"

    CombineFasta(
        input="",
        output=output,
        input_files=[missing1, missing2],
        virus_list=["VirusA", "VirusB"],
        refid_list=["RefA", "RefB"],
    ).run()

    assert output.exists()
    assert output.read_text(encoding="utf-8") == ""


@pytest.mark.xfail(reason="Mismatched metadata lengths raise KeyError instead of clear validation error", strict=False)
def test_run_rejects_mismatched_metadata_lengths_with_clear_error(tmp_path: Path) -> None:
    # Intended behavior: validate list lengths before processing and fail with ValueError.
    input1 = tmp_path / "sample.fasta"
    input1.write_text(">sampleX\nACGT\n", encoding="utf-8")
    output = tmp_path / "combined.fasta"

    combiner = CombineFasta(
        input="",
        output=output,
        input_files=[input1],
        virus_list=[],
        refid_list=[],
    )

    with pytest.raises(ValueError, match="length|metadata|virus|refid"):
        combiner.run()
