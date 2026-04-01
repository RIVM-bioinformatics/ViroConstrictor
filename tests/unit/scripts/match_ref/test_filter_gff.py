"""Unit tests for filter_gff workflow script.

This module tests the FilterGff script, which filters GFF feature records based on
reference data and remaps seqid (sequence ID) fields. Includes mock AminoExtract
stubs for hermetic testing. Tests cover happy paths, missing data, invalid types,
and filtering edge cases.
"""

import sys
import types
from argparse import ArgumentParser
from pathlib import Path

import pandas as pd
import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(PROJECT_ROOT.joinpath("ViroConstrictor/workflow")))

# Provide a minimal AminoExtract stub for hermetic unit tests.
amin_module = types.ModuleType("AminoExtract")
amin_logging = types.ModuleType("AminoExtract.logging")
amin_logging.log = object()
amin_reader = types.ModuleType("AminoExtract.reader")


class _ImportTimeSequenceReader:
    """Mock SequenceReader for AminoExtract (import-time stub)."""

    def __init__(self, logger: object, verbose: bool) -> None:
        """Initialize the mock reader.

        Parameters
        ----------
        logger : object
            Logger instance (unused in stub).
        verbose : bool
            Verbosity flag (unused in stub).
        """
        self.logger = logger
        self.verbose = verbose

    def read_gff(self, file: Path | str) -> object:
        """Placeholder for GFF reading.

        Parameters
        ----------
        file : Path | str
            GFF file path (not used in stub).

        Raises
        ------
        NotImplementedError
            Always raised in stub implementation.
        """
        raise NotImplementedError


amin_reader.SequenceReader = _ImportTimeSequenceReader
sys.modules.setdefault("AminoExtract", amin_module)
sys.modules.setdefault("AminoExtract.logging", amin_logging)
sys.modules.setdefault("AminoExtract.reader", amin_reader)

from ViroConstrictor.workflow.match_ref.scripts import filter_gff as filter_gff_module  # noqa: E402, isort:skip
from ViroConstrictor.workflow.match_ref.scripts.filter_gff import FilterGff  # noqa: E402, isort:skip


class DummyGff:
    """Mock GFF object for testing FilterGff.

    Simulates AminoExtract.reader.SequenceReader GFF object with DataFrame
    and export functionality.
    """

    def __init__(self, df: pd.DataFrame) -> None:
        """Initialize with a DataFrame.

        Parameters
        ----------
        df : pd.DataFrame
            GFF records as a DataFrame.
        """
        self.df = df
        self.exported_path: Path | None = None

    def export_gff_to_file(self, output: Path | str) -> None:
        """Export GFF to tab-separated file.

        Parameters
        ----------
        output : Path | str
            Output file path.
        """
        self.exported_path = Path(output)
        self.df.to_csv(output, sep="\t", index=False)


class FakeSequenceReader:
    """Mock SequenceReader for testing FilterGff.

    Allows injection of test GFF data without external dependencies.
    """

    def __init__(self, logger: object, verbose: bool) -> None:
        """Initialize the fake reader.

        Parameters
        ----------
        logger : object
            Logger instance (unused).
        verbose : bool
            Verbosity flag (unused).
        """
        self.logger = logger
        self.verbose = verbose
        self._gff: DummyGff | None = None

    def set_gff(self, gff: DummyGff) -> None:
        """Set the GFF data to return from read_gff.

        Parameters
        ----------
        gff : DummyGff
            Mock GFF object.
        """
        self._gff = gff

    def read_gff(self, file: Path | str) -> DummyGff:
        """Read and return the injected GFF object.

        Parameters
        ----------
        file : Path | str
            File path (unused; GFF data is pre-set).

        Returns
        -------
        DummyGff
            The injected GFF object.
        """
        assert self._gff is not None
        return self._gff


def test_add_arguments_parses_refdata_and_updatedstats() -> None:
    """Verify FilterGff.add_arguments parses refdata and updatedstats flags."""
    parser = ArgumentParser()
    FilterGff.add_arguments(parser)

    args = parser.parse_args(["--input", "in.gff", "--output", "out.gff", "--refdata", "refs.csv", "--updatedstats", "stats.csv"])
    assert args.refdata == "refs.csv"
    assert args.updatedstats == "stats.csv"


def test_run_filters_and_remaps_gff_and_updates_stats(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Verify FilterGff filters GFF records by reference mapping and updates stats.

    Tests that GFF records matching refdata are retained with remapped seqid,
    and statistics CSV is updated with output GFF file path.
    """
    refdata_csv = tmp_path / "refdata.csv"
    output_gff = tmp_path / "filtered.gff"
    updated_stats = tmp_path / "updated_stats.csv"

    pd.DataFrame(
        [
            {"seqrecord_name": "oldA", "Reference": "newA", "sample": "s1"},
            {"seqrecord_name": "oldC", "Reference": "newC", "sample": "s1"},
        ]
    ).to_csv(refdata_csv, index=False)

    source_gff = DummyGff(
        pd.DataFrame(
            [
                {"seqid": "oldA", "start": 1, "end": 10, "type": "gene"},
                {"seqid": "oldB", "start": 2, "end": 20, "type": "CDS"},
                {"seqid": "oldC", "start": 3, "end": 30, "type": "gene"},
            ]
        )
    )
    fake_reader = FakeSequenceReader(logger=object(), verbose=False)
    fake_reader.set_gff(source_gff)

    monkeypatch.setattr(filter_gff_module, "SequenceReader", lambda logger, verbose: fake_reader)

    FilterGff(
        input=tmp_path / "input.gff",
        refdata=refdata_csv,
        output=output_gff,
        updatedstats=updated_stats,
    ).run()

    out_df = pd.read_csv(output_gff, sep="\t")
    assert set(out_df["seqid"]) == {"newA", "newC"}
    assert set(out_df["type"]) == {"gene"}

    stats_df = pd.read_csv(updated_stats)
    assert "Feat_file" in stats_df.columns
    assert set(stats_df["Feat_file"]) == {str(output_gff.resolve())}


def test_run_with_missing_ref_columns_raises_key_error(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Verify that missing expected columns in refdata CSV raises KeyError."""
    refdata_csv = tmp_path / "refdata_invalid.csv"
    pd.DataFrame([{"wrong": "oldA", "also_wrong": "newA"}]).to_csv(refdata_csv, index=False)

    fake_reader = FakeSequenceReader(logger=object(), verbose=False)
    fake_reader.set_gff(DummyGff(pd.DataFrame([{"seqid": "oldA", "start": 1, "end": 10}])))
    monkeypatch.setattr(filter_gff_module, "SequenceReader", lambda logger, verbose: fake_reader)

    with pytest.raises(KeyError):
        FilterGff(
            input=tmp_path / "input.gff",
            refdata=refdata_csv,
            output=tmp_path / "filtered.gff",
            updatedstats=tmp_path / "updated.csv",
        ).run()


def test_run_rejects_invalid_input_type(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Verify that invalid input type (not Path/str) raises AssertionError."""
    refdata_csv = tmp_path / "refdata.csv"
    pd.DataFrame([{"seqrecord_name": "oldA", "Reference": "newA"}]).to_csv(refdata_csv, index=False)

    fake_reader = FakeSequenceReader(logger=object(), verbose=False)
    fake_reader.set_gff(DummyGff(pd.DataFrame([{"seqid": "oldA", "start": 1, "end": 10}])))
    monkeypatch.setattr(filter_gff_module, "SequenceReader", lambda logger, verbose: fake_reader)

    script = FilterGff(
        input=tmp_path / "input.gff",
        refdata=refdata_csv,
        output=tmp_path / "filtered.gff",
        updatedstats=tmp_path / "updated.csv",
    )
    script.input = ["not", "a", "path"]  # type: ignore[assignment]

    with pytest.raises(AssertionError, match="Input"):
        script.run()


@pytest.mark.xfail(
    strict=True,
    reason="Intended behavior: no matching GFF records should raise a clear error rather than writing empty output.",
)
def test_run_raises_when_no_gff_records_match_reference_data(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Verify that no matching GFF records raises ValueError with clear message."""
    refdata_csv = tmp_path / "refdata.csv"
    pd.DataFrame([{"seqrecord_name": "oldA", "Reference": "newA"}]).to_csv(refdata_csv, index=False)

    fake_reader = FakeSequenceReader(logger=object(), verbose=False)
    fake_reader.set_gff(DummyGff(pd.DataFrame([{"seqid": "unmatched", "start": 1, "end": 10}])))
    monkeypatch.setattr(filter_gff_module, "SequenceReader", lambda logger, verbose: fake_reader)

    with pytest.raises(ValueError, match="No matching GFF entries"):
        FilterGff(
            input=tmp_path / "input.gff",
            refdata=refdata_csv,
            output=tmp_path / "filtered.gff",
            updatedstats=tmp_path / "updated.csv",
        ).run()
