"""Unit tests for the `GroupRefs` workflow script.

These tests exercise the grouping logic used by the match-ref workflow's
`group_refs` script. Tests verify segmented and non-segmented mode, argument
parsing, and comprehensive error handling for invalid inputs.
"""

import runpy
import sys
from argparse import ArgumentParser
from pathlib import Path
from typing import List, cast

import pandas as pd
import pytest
from Bio import SeqIO

PROJECT_ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(PROJECT_ROOT.joinpath("ViroConstrictor/workflow")))
from ViroConstrictor.workflow.match_ref.scripts.group_refs import GroupRefs  # noqa: E402, isort:skip


def test_group_refs_segmented_mode(tmp_path: Path) -> None:
    """
    Verify grouping when multiple segmented reference files are provided.

    Parameters
    ----------
    tmp_path : Path
        Temporary directory provided by pytest.

    Returns
    -------
    None
    """

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
    """
    Ensure a single non-segmented reference keeps its original sequence id.

    Parameters
    ----------
    tmp_path : Path
        Temporary directory provided by pytest.

    Returns
    -------
    None
    """

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
    """
    Passing a non-list for `input_refs` should be rejected.

    Parameters
    ----------
    tmp_path : Path
        Temporary directory provided by pytest.

    Returns
    -------
    None
    """

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

    with pytest.raises((AssertionError, TypeError, ValueError)):
        grouper.run()


def test_group_refs_invalid_input_stats_type_raises(tmp_path: Path) -> None:
    """Passing a non-list for `input_stats` should be rejected.

    Parameters
    ----------
    tmp_path : Path
        Temporary directory provided by pytest.

    Returns
    -------
    None
    """

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

    with pytest.raises((AssertionError, TypeError, ValueError)):
        grouper.run()


def test_group_refs_missing_input_ref_file_raises(tmp_path: Path) -> None:
    """
    A missing reference file should raise FileNotFoundError.

    Parameters
    ----------
    tmp_path : Path
        Temporary directory provided by pytest.

    Returns
    -------
    None
    """

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
    """A missing stats CSV file should raise FileNotFoundError.

    Parameters
    ----------
    tmp_path : Path
        Temporary directory provided by pytest.

    Returns
    -------
    None
    """

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
    """If the stats CSV lacks a `Reference` column, processing should fail.

    Parameters
    ----------
    tmp_path : Path
        Temporary directory provided by pytest.

    Returns
    -------
    None
    """

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

    with pytest.raises((KeyError, ValueError), match="Reference"):
        grouper.run()


@pytest.mark.xfail(
    strict=True,
    reason="Current implementation raises IndexError; intended behavior is a clear validation error for empty FASTA input.",
)
def test_group_refs_empty_fasta_without_records_raises_value_error(tmp_path: Path) -> None:
    """An empty FASTA file (no records) should be rejected with a validation error.

    Parameters
    ----------
    tmp_path : Path
        Temporary directory provided by pytest.

    Returns
    -------
    None
    """

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

    with pytest.raises(ValueError, match="fasta|reference|record"):
        grouper.run()


def test_group_refs_add_arguments_requires_expected_flags() -> None:
    """Ensure parser wiring exposes both base and script-specific required flags."""

    parser = ArgumentParser()
    GroupRefs.add_arguments(parser)

    args = parser.parse_args(
        [
            "--input",
            "unused",
            "--output",
            "out.fasta",
            "--input_refs",
            "in1.fasta",
            "in2.fasta",
            "--input_stats",
            "in1.csv",
            "in2.csv",
            "--output_stats",
            "out.csv",
            "--sample",
            "sample-a",
        ]
    )

    assert args.input == "unused"
    assert args.output == "out.fasta"
    assert args.input_refs == ["in1.fasta", "in2.fasta"]
    assert args.input_stats == ["in1.csv", "in2.csv"]
    assert args.output_stats == "out.csv"
    assert args.sample == "sample-a"


def test_group_refs_main_block_invokes_class_main(monkeypatch: pytest.MonkeyPatch) -> None:
    """Verify running the module as a script triggers GroupRefs.main via the main guard."""

    import helpers.base_script_class as base_script_class

    called: dict[str, str] = {}

    def fake_main(cls: type[object]) -> None:
        called["class_name"] = cls.__name__

    monkeypatch.setattr(base_script_class.BaseScript, "main", classmethod(fake_main))

    runpy.run_path(
        str(PROJECT_ROOT / "ViroConstrictor" / "workflow" / "match_ref" / "scripts" / "group_refs.py"),
        run_name="__main__",
    )

    assert called == {"class_name": "GroupRefs"}


@pytest.mark.xfail(
    strict=True,
    reason="Current implementation silently writes rows with NaN Reference instead of rejecting unmatched stats rows.",
)
def test_group_refs_rejects_stats_rows_that_do_not_match_any_reference(tmp_path: Path) -> None:
    """Intended behavior: fail when no stats row can be mapped to grouped reference records."""

    input_ref = str(PROJECT_ROOT / "tests" / "unit" / "data" / "group_refs" / "HA_INF1024_B34_4312201729_squiggle_best_ref.fasta")
    stats_file = tmp_path / "unmatched_stats.csv"
    stats_file.write_text("Reference,Mapped Reads\ntotally_unrelated_ref,123\n", encoding="utf-8")

    grouper = GroupRefs(
        input="unused",
        input_refs=cast(List[Path | str], [input_ref]),
        input_stats=cast(List[Path | str], [stats_file]),
        output=tmp_path / "grouped_refs.fasta",
        output_stats=tmp_path / "grouped_stats.csv",
        sample="sample",
    )

    with pytest.raises(ValueError, match="Reference"):
        grouper.run()


@pytest.mark.xfail(
    strict=True,
    reason="Current implementation raises IndexError for empty input_refs; intended behavior is a clear validation error.",
)
def test_group_refs_empty_input_refs_list_raises_value_error(tmp_path: Path) -> None:
    """Intended behavior: reject an empty list of reference files with ValueError."""

    input_stat = str(PROJECT_ROOT / "tests" / "unit" / "data" / "group_refs" / "HA_INF1024_B34_4312201729_squiggle_best_ref.csv")

    grouper = GroupRefs(
        input="unused",
        input_refs=cast(List[Path | str], []),
        input_stats=cast(List[Path | str], [input_stat]),
        output=tmp_path / "grouped_refs.fasta",
        output_stats=tmp_path / "grouped_stats.csv",
        sample="sample",
    )

    with pytest.raises(ValueError, match="reference"):
        grouper.run()
