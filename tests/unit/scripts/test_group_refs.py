import sys
from pathlib import Path
from typing import List, cast

import pandas as pd
import pytest
from Bio import SeqIO

PROJECT_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PROJECT_ROOT.joinpath("ViroConstrictor/workflow")))
from ViroConstrictor.workflow.match_ref.scripts.group_refs import GroupRefs  # isort:skip


def test_group_refs_segmented_mode(tmp_path: Path) -> None:
    _input_refs = [
        str(PROJECT_ROOT / "tests" / "unit" / "data" / "group_refs" / "HA_INF1024_B34_4312201729_squiggle_best_ref.fasta"),
        str(PROJECT_ROOT / "tests" / "unit" / "data" / "group_refs" / "MP_INF1024_B34_4312201729_squiggle_best_ref.fasta"),
    ]
    _input_stats = [
        str(PROJECT_ROOT / "tests" / "unit" / "data" / "group_refs" / "HA_INF1024_B34_4312201729_squiggle_best_ref.csv"),
        str(PROJECT_ROOT / "tests" / "unit" / "data" / "group_refs" / "MP_INF1024_B34_4312201729_squiggle_best_ref.csv"),
    ]

    output = tmp_path / "output_refs.fasta"
    output_stats = tmp_path / "output_stats.csv"
    sample_name = "INF1024_B34_4312201729_squiggle"

    grouper = GroupRefs(
        input="empty",
        input_refs=cast(List[Path | str], _input_refs),
        input_stats=cast(List[Path | str], _input_stats),
        output=output,
        output_stats=output_stats,
        sample=sample_name,
    )

    grouper.run()

    assert output.exists()
    assert output_stats.exists()

    records = list(SeqIO.parse(output, "fasta"))
    assert len(records) == 2
    assert {record.id for record in records} == {"HA", "MP"}

    stats = pd.read_csv(output_stats)
    assert len(stats) == 2
    assert set(stats["Reference"]) == {"HA", "MP"}
    assert set(stats["sample"]) == {sample_name}
    assert set(stats["seqrecord_id"]) == {"HA", "MP"}
    assert set(stats["Reference_file"]) == {str(output.resolve())}


def test_group_refs_non_segmented_mode_keeps_id(tmp_path: Path) -> None:
    input_ref = str(PROJECT_ROOT / "tests" / "unit" / "data" / "group_refs" / "HA_INF1024_B34_4312201729_squiggle_best_ref.fasta")
    input_stat = str(PROJECT_ROOT / "tests" / "unit" / "data" / "group_refs" / "HA_INF1024_B34_4312201729_squiggle_best_ref.csv")

    output = tmp_path / "output_refs_single.fasta"
    output_stats = tmp_path / "output_stats_single.csv"

    grouper = GroupRefs(
        input="empty",
        input_refs=cast(List[Path | str], [input_ref]),
        input_stats=cast(List[Path | str], [input_stat]),
        output=output,
        output_stats=output_stats,
        sample="single_sample",
    )

    grouper.run()

    records = list(SeqIO.parse(output, "fasta"))
    assert len(records) == 1
    assert records[0].id == "A.HA_10"

    stats = pd.read_csv(output_stats)
    assert len(stats) == 1
    assert stats.loc[0, "Reference"] == "A.HA_10"
    assert stats.loc[0, "sample"] == "single_sample"


def test_group_refs_invalid_input_refs_type_raises(tmp_path: Path) -> None:
    output = tmp_path / "output_refs.fasta"
    output_stats = tmp_path / "output_stats.csv"
    input_stat = str(PROJECT_ROOT / "tests" / "unit" / "data" / "group_refs" / "HA_INF1024_B34_4312201729_squiggle_best_ref.csv")

    grouper = GroupRefs(
        input="empty",
        input_refs=cast(List[Path | str], "not-a-list"),
        input_stats=cast(List[Path | str], [input_stat]),
        output=output,
        output_stats=output_stats,
        sample="sample",
    )

    with pytest.raises(AssertionError, match="Input_refs should be a list"):
        grouper.run()


def test_group_refs_invalid_input_stats_type_raises(tmp_path: Path) -> None:
    output = tmp_path / "output_refs.fasta"
    output_stats = tmp_path / "output_stats.csv"
    input_ref = str(PROJECT_ROOT / "tests" / "unit" / "data" / "group_refs" / "HA_INF1024_B34_4312201729_squiggle_best_ref.fasta")

    grouper = GroupRefs(
        input="empty",
        input_refs=cast(List[Path | str], [input_ref]),
        input_stats=cast(List[Path | str], "not-a-list"),
        output=output,
        output_stats=output_stats,
        sample="sample",
    )

    with pytest.raises(AssertionError, match="Input_stats should be a list"):
        grouper.run()


def test_group_refs_missing_input_ref_file_raises(tmp_path: Path) -> None:
    output = tmp_path / "output_refs.fasta"
    output_stats = tmp_path / "output_stats.csv"
    missing_ref = tmp_path / "missing_ref.fasta"
    input_stat = str(PROJECT_ROOT / "tests" / "unit" / "data" / "group_refs" / "HA_INF1024_B34_4312201729_squiggle_best_ref.csv")

    grouper = GroupRefs(
        input="empty",
        input_refs=cast(List[Path | str], [str(missing_ref)]),
        input_stats=cast(List[Path | str], [input_stat]),
        output=output,
        output_stats=output_stats,
        sample="sample",
    )

    with pytest.raises(FileNotFoundError):
        grouper.run()


def test_group_refs_missing_input_stats_file_raises(tmp_path: Path) -> None:
    output = tmp_path / "output_refs.fasta"
    output_stats = tmp_path / "output_stats.csv"
    missing_stats = tmp_path / "missing_stats.csv"
    input_ref = str(PROJECT_ROOT / "tests" / "unit" / "data" / "group_refs" / "HA_INF1024_B34_4312201729_squiggle_best_ref.fasta")

    grouper = GroupRefs(
        input="empty",
        input_refs=cast(List[Path | str], [input_ref]),
        input_stats=cast(List[Path | str], [str(missing_stats)]),
        output=output,
        output_stats=output_stats,
        sample="sample",
    )

    with pytest.raises(FileNotFoundError):
        grouper.run()


def test_group_refs_stats_without_reference_column_raises(tmp_path: Path) -> None:
    input_ref = str(PROJECT_ROOT / "tests" / "unit" / "data" / "group_refs" / "HA_INF1024_B34_4312201729_squiggle_best_ref.fasta")
    invalid_stats = tmp_path / "invalid_stats.csv"
    output = tmp_path / "output_refs.fasta"
    output_stats = tmp_path / "output_stats.csv"

    pd.DataFrame(
        {
            "Mapped Reads": [1],
            "Avg. Mismatches per Read": [0.1],
            "Avg. Sequence Identity": [0.9],
        }
    ).to_csv(invalid_stats, index=False)

    grouper = GroupRefs(
        input="empty",
        input_refs=cast(List[Path | str], [input_ref]),
        input_stats=cast(List[Path | str], [invalid_stats]),
        output=output,
        output_stats=output_stats,
        sample="sample",
    )

    with pytest.raises(KeyError, match="Reference"):
        grouper.run()


def test_group_refs_empty_fasta_without_records_raises(tmp_path: Path) -> None:
    empty_ref = tmp_path / "empty.fasta"
    input_stat = str(PROJECT_ROOT / "tests" / "unit" / "data" / "group_refs" / "HA_INF1024_B34_4312201729_squiggle_best_ref.csv")
    output = tmp_path / "output_refs.fasta"
    output_stats = tmp_path / "output_stats.csv"

    empty_ref.write_text("", encoding="utf-8")

    grouper = GroupRefs(
        input="empty",
        input_refs=cast(List[Path | str], [str(empty_ref)]),
        input_stats=cast(List[Path | str], [input_stat]),
        output=output,
        output_stats=output_stats,
        sample="sample",
    )

    with pytest.raises(IndexError):
        grouper.run()
