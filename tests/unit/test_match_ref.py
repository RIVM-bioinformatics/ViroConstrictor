"""Tests for :mod:`ViroConstrictor.match_ref`.

This module validates the dataframe merge helpers and the high-level match-ref
workflow orchestration, including dry-run handling, failure reporting, and
sample table updates.
"""

from types import SimpleNamespace
from typing import Any, cast

import pandas as pd
import pytest

import ViroConstrictor.match_ref as match_ref


def _make_parsed_inputs(samples_df: pd.DataFrame) -> SimpleNamespace:
    """Build a minimal parsed-input namespace for match-ref tests.

    Parameters
    ----------
    samples_df : pandas.DataFrame
        Sample table to attach to the namespace.

    Returns
    -------
    types.SimpleNamespace
        Namespace with the attributes required by ``process_match_ref``.
    """
    return SimpleNamespace(
        workdir="/tmp/work",
        input_path="/tmp/input",
        exec_start_path="/tmp/start",
        user_config={"x": "y"},
        samples_df=samples_df,
        samples_dict={},
    )


def test_replace_sets_to_singular_values_replaces_target_columns() -> None:
    """Test replace_sets_to_singular_values extracts single items from set-valued columns.
    
    Verifies that set columns are converted to their single string values for specified
    column names.
    """
    df = pd.DataFrame(
        {
            "sample": ["S1", "S2"],
            "Reference_file": [{"ref_s1"}, {"ref_s2"}],
            "Primer_file": [{"primer_s1"}, {"primer_s2"}],
            "Feat_file": [{"feat_s1"}, {"feat_s2"}],
        }
    )

    out = match_ref.replace_sets_to_singular_values(df, ["Reference_file", "Primer_file", "Feat_file"])

    assert out.loc[0, "Reference_file"] == "ref_s1"
    assert out.loc[1, "Primer_file"] == "primer_s2"
    assert out.loc[0, "Feat_file"] == "feat_s1"


def test_replacement_merge_dataframe_on_cols_respects_none_and_missing_sample() -> None:
    """Test replacement_merge_dataframe_on_cols merges override values while respecting NONE.
    
    Verifies that:
    - New values from override_df replace old values in original_df by sample.
    - "NONE" values are preserved and not replaced.
    - Missing samples in override_df retain their original values.
    """
    original_df = pd.DataFrame(
        {
            "SAMPLE": [1, 2, 3],
            "REFERENCE": ["old_ref_1", "NONE", "old_ref_3"],
            "PRIMERS": ["old_primer_1", "old_primer_2", "NONE"],
            "FEATURES": ["old_feat_1", "old_feat_2", "old_feat_3"],
        }
    )
    override_df = pd.DataFrame(
        {
            "sample": ["1", "2"],
            "Reference_file": ["new_ref_1", "new_ref_2"],
            "Primer_file": ["new_primer_1", "new_primer_2"],
            "Feat_file": ["new_feat_1", "new_feat_2"],
        }
    )

    out = match_ref.replacement_merge_dataframe_on_cols(
        original_df,
        override_df,
        ["REFERENCE", "PRIMERS", "FEATURES"],
        ["Reference_file", "Primer_file", "Feat_file"],
    )

    assert out.loc[0, "REFERENCE"] == "new_ref_1"
    assert out.loc[0, "PRIMERS"] == "new_primer_1"
    assert out.loc[1, "REFERENCE"] == "NONE"
    assert out.loc[1, "FEATURES"] == "new_feat_2"
    assert out.loc[2, "REFERENCE"] == "old_ref_3"
    assert out.loc[2, "PRIMERS"] == "NONE"


def test_process_match_ref_failure_writes_report_and_exits(monkeypatch: pytest.MonkeyPatch) -> None:
    """Test process_match_ref writes failure report and exits with code 1 on workflow failure.
    
    Verifies that:
    - Match-ref workflow stage ("MR") is invoked.
    - On workflow failure (return False), a failure report is written.
    - Process exits with code 1.
    """
    parsed_inputs = _make_parsed_inputs(
        pd.DataFrame(
            {
                "SAMPLE": ["S1"],
                "REFERENCE": ["old_ref"],
                "PRIMERS": ["old_primer"],
                "FEATURES": ["old_feat"],
            }
        )
    )
    used_workflow_config = SimpleNamespace(
        resource_settings=SimpleNamespace(),
        output_settings=SimpleNamespace(dryrun=False),
    )

    called: dict[str, Any] = {}

    def fake_run_snakemake_workflow(*args: Any, **kwargs: Any) -> tuple[bool, SimpleNamespace]:
        """Mock workflow executor that records stage/scheduler and returns failure status."""
        called["stage"] = kwargs.get("stage")
        called["scheduler"] = kwargs.get("scheduler")
        return False, used_workflow_config

    def fake_write_report(*args: Any, **kwargs: Any) -> None:
        """Mock report writer that records all report invocation arguments."""
        called["report_args"] = args
        called["report_kwargs"] = kwargs

    monkeypatch.setattr(match_ref, "run_snakemake_workflow", fake_run_snakemake_workflow)
    monkeypatch.setattr(match_ref, "WriteReport", fake_write_report)

    with pytest.raises(SystemExit) as exc_info:
        match_ref.process_match_ref(
            parsed_inputs=cast(match_ref.CLIparser, parsed_inputs),
            scheduler=cast(match_ref.Scheduler, object()),
        )

    assert exc_info.value.code == 1
    assert called["stage"] == "MR"
    assert called["report_args"][-1] == "Failed"


def test_process_match_ref_dryrun_returns_unmodified_input(monkeypatch: pytest.MonkeyPatch) -> None:
    """Test process_match_ref skips pickle loading and returns input unchanged in dryrun mode.
    
    Verifies that when dryrun is True:
    - The workflow is run but result pickle is not loaded.
    - The original parsed_inputs object is returned unmodified.
    - The samples dataframe retains its original columns.
    """
    parsed_inputs = _make_parsed_inputs(
        pd.DataFrame(
            {
                "SAMPLE": ["S1"],
                "REFERENCE": ["old_ref"],
                "PRIMERS": ["old_primer"],
                "FEATURES": ["old_feat"],
            }
        )
    )

    used_workflow_config = SimpleNamespace(
        resource_settings=SimpleNamespace(),
        output_settings=SimpleNamespace(dryrun=True),
    )

    monkeypatch.setattr(match_ref, "run_snakemake_workflow", lambda *args, **kwargs: (True, used_workflow_config))

    def _should_not_read_pickle(_path: str) -> pd.DataFrame:
        """Sentinel function that raises if called, used to verify dryrun behavior."""
        raise AssertionError("read_pickle should not be called in dryrun mode")

    monkeypatch.setattr(match_ref.pd, "read_pickle", _should_not_read_pickle)

    out = match_ref.process_match_ref(
        parsed_inputs=cast(match_ref.CLIparser, parsed_inputs),
        scheduler=cast(match_ref.Scheduler, object()),
    )

    assert out is parsed_inputs
    assert list(out.samples_df.columns) == ["SAMPLE", "REFERENCE", "PRIMERS", "FEATURES"]


def test_process_match_ref_success_updates_samples_and_sets_index(monkeypatch: pytest.MonkeyPatch) -> None:
    """Test process_match_ref updates samples from match-ref results and sets SAMPLE as index.
    
    Verifies that on successful workflow execution:
    - Sample columns (REFERENCE, PRIMERS, FEATURES) are updated from pickle results.
    - The SAMPLE column is set as the dataframe index.
    - The samples_dict is populated with sample-level reference mappings.
    """
    parsed_inputs = _make_parsed_inputs(
        pd.DataFrame(
            {
                "SAMPLE": ["S1", "S2"],
                "REFERENCE": ["old_ref_1", "NONE"],
                "PRIMERS": ["old_primer_1", "NONE"],
                "FEATURES": ["old_feat_1", "old_feat_2"],
            }
        )
    )
    used_workflow_config = SimpleNamespace(
        resource_settings=SimpleNamespace(),
        output_settings=SimpleNamespace(dryrun=False),
    )

    result_df = pd.DataFrame(
        {
            "sample": ["S1", "S1", "S2"],
            "Reference": ["R1", "R1", "R2"],
            "Reference_file": ["new_ref_1", "new_ref_1", "new_ref_2"],
        }
    )

    monkeypatch.setattr(match_ref, "run_snakemake_workflow", lambda *args, **kwargs: (True, used_workflow_config))
    monkeypatch.setattr(match_ref.pd, "read_pickle", lambda _path: result_df)

    out = match_ref.process_match_ref(
        parsed_inputs=cast(match_ref.CLIparser, parsed_inputs),
        scheduler=cast(match_ref.Scheduler, object()),
    )

    assert out.samples_df.index.name == "SAMPLE"
    assert out.samples_df.loc["S1", "REFERENCE"] == "new_ref_1"
    assert out.samples_df.loc["S2", "REFERENCE"] == "NONE"
    assert out.samples_df.loc["S1", "PRIMERS"] == "NONE"
    assert out.samples_df.loc["S1", "FEATURES"] == "NONE"
    assert out.samples_dict["S1"]["REFERENCE"] == "new_ref_1"


def test_process_match_ref_fills_nan_before_grouping(monkeypatch: pytest.MonkeyPatch) -> None:
    """Test process_match_ref converts NaN and None to "NONE" before updating samples.
    
    Verifies that missing values (None, pd.NA) in match-ref results are replaced with
    the "NONE" sentinel string in the final samples dataframe.
    """
    parsed_inputs = _make_parsed_inputs(
        pd.DataFrame(
            {
                "SAMPLE": ["S1"],
                "REFERENCE": ["old_ref_1"],
                "PRIMERS": ["old_primer_1"],
                "FEATURES": ["old_feat_1"],
            }
        )
    )
    used_workflow_config = SimpleNamespace(
        resource_settings=SimpleNamespace(),
        output_settings=SimpleNamespace(dryrun=False),
    )

    result_df = pd.DataFrame(
        {
            "sample": ["S1"],
            "Reference": ["R1"],
            "Reference_file": ["new_ref_1"],
            "Primer_file": [None],
            "Feat_file": [pd.NA],
        }
    )

    monkeypatch.setattr(match_ref, "run_snakemake_workflow", lambda *args, **kwargs: (True, used_workflow_config))
    monkeypatch.setattr(match_ref.pd, "read_pickle", lambda _path: result_df)

    out = match_ref.process_match_ref(
        parsed_inputs=cast(match_ref.CLIparser, parsed_inputs),
        scheduler=cast(match_ref.Scheduler, object()),
    )

    assert out.samples_df.loc["S1", "PRIMERS"] == "NONE"
    assert out.samples_df.loc["S1", "FEATURES"] == "NONE"
