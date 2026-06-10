"""Unit tests for filter_bed workflow script.

This module tests the FilterBed script, which filters a BED file (primer coordinates)
based on reference mappings, remaps reference IDs, and updates statistics. Tests cover
happy paths, edge cases (empty matches, missing files), and error handling.
"""

import sys
from argparse import ArgumentParser
from pathlib import Path

import pandas as pd
import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(PROJECT_ROOT.joinpath("ViroConstrictor/workflow")))
from ViroConstrictor.workflow.match_ref.scripts.filter_bed import FilterBed  # noqa: E402, isort:skip


def test_add_arguments_parses_refdata_and_updatedstats() -> None:
    """Verify that FilterBed.add_arguments parses refdata and updatedstats flags."""
    parser = ArgumentParser()
    FilterBed.add_arguments(parser)

    args = parser.parse_args(["--input", "in.bed", "--output", "out.bed", "--refdata", "refs.csv", "--updatedstats", "stats.csv"])
    assert args.refdata == "refs.csv"
    assert args.updatedstats == "stats.csv"


def test_run_filters_and_maps_bed_and_updates_stats(tmp_path: Path) -> None:
    """Verify FilterBed filters BED records by reference mapping and updates stats.

    Tests that BED records matching refdata are retained with remapped reference IDs,
    and statistics CSV is updated with the output BED file path.
    """
    input_bed = tmp_path / "primers.bed"
    refdata_csv = tmp_path / "refdata.csv"
    output_bed = tmp_path / "filtered.bed"
    updated_stats = tmp_path / "updated_stats.csv"

    input_bed.write_text(
        "oldA\t1\t10\tamp1\t0\t+\noldB\t2\t20\tamp2\t0\t-\noldC\t3\t30\tamp3\t0\t+\n",
        encoding="utf-8",
    )
    pd.DataFrame(
        [
            {"seqrecord_name": "oldA", "Reference": "newA", "sample": "s1"},
            {"seqrecord_name": "oldC", "Reference": "newC", "sample": "s1"},
        ]
    ).to_csv(refdata_csv, index=False)

    FilterBed(input=input_bed, refdata=refdata_csv, output=output_bed, updatedstats=updated_stats).run()

    out_df = pd.read_csv(output_bed, sep="\t", header=None)
    assert out_df.shape == (2, 6)
    assert set(out_df.iloc[:, 0]) == {"newA", "newC"}
    assert set(out_df.iloc[:, 3]) == {"amp1", "amp3"}

    stats_df = pd.read_csv(updated_stats)
    assert "Primer_file" in stats_df.columns
    assert set(stats_df["Primer_file"]) == {str(output_bed.resolve())}


def test_run_when_no_bed_rows_match_writes_empty_output(tmp_path: Path) -> None:
    """Verify FilterBed writes empty BED file when no records match reference data."""
    input_bed = tmp_path / "primers.bed"
    refdata_csv = tmp_path / "refdata.csv"
    output_bed = tmp_path / "filtered.bed"
    updated_stats = tmp_path / "updated_stats.csv"

    input_bed.write_text("oldB\t2\t20\tamp2\t0\t-\n", encoding="utf-8")
    pd.DataFrame([{"seqrecord_name": "oldA", "Reference": "newA"}]).to_csv(refdata_csv, index=False)

    FilterBed(input=input_bed, refdata=refdata_csv, output=output_bed, updatedstats=updated_stats).run()

    assert output_bed.exists()
    assert output_bed.read_text(encoding="utf-8") == ""

    stats_df = pd.read_csv(updated_stats)
    assert len(stats_df) == 1


def test_run_with_missing_ref_columns_raises_key_error(tmp_path: Path) -> None:
    """Verify that missing expected columns in refdata CSV raises KeyError."""
    input_bed = tmp_path / "primers.bed"
    refdata_csv = tmp_path / "refdata.csv"

    input_bed.write_text("oldA\t1\t10\tamp1\t0\t+\n", encoding="utf-8")
    pd.DataFrame([{"wrong": "oldA", "also_wrong": "newA"}]).to_csv(refdata_csv, index=False)

    with pytest.raises(KeyError):
        FilterBed(
            input=input_bed,
            refdata=refdata_csv,
            output=tmp_path / "filtered.bed",
            updatedstats=tmp_path / "updated.csv",
        ).run()


def test_run_with_missing_bed_file_raises_file_not_found(tmp_path: Path) -> None:
    """Verify that a missing input BED file raises FileNotFoundError."""
    refdata_csv = tmp_path / "refdata.csv"
    pd.DataFrame([{"seqrecord_name": "oldA", "Reference": "newA"}]).to_csv(refdata_csv, index=False)

    with pytest.raises(FileNotFoundError):
        FilterBed(
            input=tmp_path / "does_not_exist.bed",
            refdata=refdata_csv,
            output=tmp_path / "filtered.bed",
            updatedstats=tmp_path / "updated.csv",
        ).run()


@pytest.mark.xfail(
    strict=True,
    reason="Intended behavior: malformed BED records should raise ValueError instead of being silently tolerated.",
)
def test_run_rejects_malformed_bed_rows(tmp_path: Path) -> None:
    """Verify that malformed BED rows (missing columns) are rejected with ValueError."""
    input_bed = tmp_path / "malformed.bed"
    refdata_csv = tmp_path / "refdata.csv"

    # Missing the strand column.
    input_bed.write_text("oldA\t1\t10\tamp1\t0\n", encoding="utf-8")
    pd.DataFrame([{"seqrecord_name": "oldA", "Reference": "newA"}]).to_csv(refdata_csv, index=False)

    with pytest.raises(ValueError, match="BED|columns|format"):
        FilterBed(
            input=input_bed,
            refdata=refdata_csv,
            output=tmp_path / "filtered.bed",
            updatedstats=tmp_path / "updated.csv",
        ).run()
