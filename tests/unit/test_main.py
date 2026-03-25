"""Test main entrypoint orchestration and preset warning handling.

These tests cover warning aggregation, preset logging, and the high-level main
control flow that coordinates parser, workflow, and reporting helpers.
"""

from __future__ import annotations

from types import SimpleNamespace
from typing import Any

import pandas as pd
import pytest

import ViroConstrictor.__main__ as vc_main


def test_get_preset_warning_list_returns_empty_lists_when_no_matches() -> None:
    """Test that get_preset_warning_list returns empty warnings when presets match perfectly.
    
    Verifies that no fallback or score-based warnings are generated when the preset
    score is high and the preset is not the DEFAULT fallback.
    """
    df = pd.DataFrame(
        {
            "VIRUS": ["SARS-CoV-2"],
            "PRESET": ["SARSCOV2"],
            "SAMPLE": ["S1"],
            "PRESET_SCORE": [0.9],
        }
    )

    fallback_warnings, score_warnings = vc_main.get_preset_warning_list(df)

    assert fallback_warnings == []
    assert score_warnings == []


def test_get_preset_warning_list_builds_score_and_fallback_warnings() -> None:
    """Test that get_preset_warning_list generates score and fallback warnings appropriately.
    
    Verifies that:
    - Samples with low preset match scores (<80%) generate score-based warnings.
    - Samples that default to the DEFAULT preset generate fallback warnings.
    """
    df = pd.DataFrame(
        {
            "VIRUS": ["SARS", "SARS", "mystery-virus"],
            "PRESET": ["SARSCOV2", "SARSCOV2", "DEFAULT"],
            "SAMPLE": ["S1", "S2", "S3"],
            "PRESET_SCORE": [0.5, 0.5, 0.0],
        }
    )

    fallback_warnings, score_warnings = vc_main.get_preset_warning_list(df)

    assert len(score_warnings) == 1
    assert "less than 80% certainty" in score_warnings[0]
    assert "S1" in score_warnings[0]
    assert "S2" in score_warnings[0]
    assert "50%" in score_warnings[0]

    assert len(fallback_warnings) == 1
    assert "could however not be used to determine a preset" in fallback_warnings[0]
    assert "DEFAULT" in fallback_warnings[0]
    assert "S3" in fallback_warnings[0]


def test_show_preset_warnings_with_per_sample_disable_column(monkeypatch: pytest.MonkeyPatch) -> None:
    """Test show_preset_warnings respects per-sample DISABLE-PRESETS column.
    
    Verifies that warnings are filtered based on individual sample disable flags
    when the DISABLE-PRESETS column is present in the samples dataframe.
    """
    samples_df = pd.DataFrame(
        {
            "DISABLE-PRESETS": [False, True],
        },
        index=["sample_a", "sample_b"],
    )

    warnings = ["warning for sample_a", "warning for sample_b"]
    fallbacks = ["fallback for sample_a", "fallback for sample_b"]
    logged_messages: list[str] = []
    monkeypatch.setattr(vc_main.log, "warning", lambda message: logged_messages.append(message))

    vc_main.show_preset_warnings(warnings, fallbacks, disabled=False, samples_df=samples_df)

    assert logged_messages == ["warning for sample_a", "fallback for sample_a"]


def test_show_preset_warnings_uses_global_disable_flag_without_column(monkeypatch: pytest.MonkeyPatch) -> None:
    """Test show_preset_warnings uses global disable flag when DISABLE-PRESETS column absent.
    
    Verifies that warnings are shown or suppressed based on the global disabled parameter
    when the samples dataframe lacks a per-sample DISABLE-PRESETS column.
    """
    samples_df = pd.DataFrame({"SAMPLE": ["S1"]})
    logged_messages: list[str] = []
    monkeypatch.setattr(vc_main.log, "warning", lambda message: logged_messages.append(message))

    vc_main.show_preset_warnings(
        warnings=["score warning"],
        fallbacks=["fallback warning"],
        disabled=False,
        samples_df=samples_df,
    )

    assert logged_messages == ["score warning", "fallback warning"]


def test_show_preset_warnings_skips_when_disabled(monkeypatch: pytest.MonkeyPatch) -> None:
    """Test show_preset_warnings suppresses all warnings when disabled is True.
    
    Verifies that both score and fallback warnings are skipped when the disabled flag is set.
    """
    samples_df = pd.DataFrame({"SAMPLE": ["S1"]})
    logged_messages: list[str] = []
    monkeypatch.setattr(vc_main.log, "warning", lambda message: logged_messages.append(message))

    vc_main.show_preset_warnings(
        warnings=["score warning"],
        fallbacks=["fallback warning"],
        disabled=True,
        samples_df=samples_df,
    )

    assert logged_messages == []


def _build_parsed_input(*, skip_updates: bool, disable_presets: bool, match_ref: bool) -> SimpleNamespace:
    """Build parsed input state for main entrypoint tests.

    Parameters
    ----------
    skip_updates : bool
        Flag indicating whether update steps should be skipped.
    disable_presets : bool
        Flag indicating whether preset warnings should be suppressed.
    match_ref : bool
        Controls whether the parsed samples dataframe requires match-ref work.

    Returns
    -------
    SimpleNamespace
        Parsed input object with the attributes accessed by main().
    """
    return SimpleNamespace(
        samples_df=pd.DataFrame({"MATCH-REF": [match_ref]}),
        flags=SimpleNamespace(skip_updates=skip_updates, disable_presets=disable_presets),
        user_config={"profile": "x"},
        scheduler="LOCAL",
        workdir="/tmp/workdir",
        input_path="/tmp/input",
        exec_start_path="/tmp/start",
    )


def test_main_happy_path_runs_update_match_ref_and_exits_zero(monkeypatch: pytest.MonkeyPatch) -> None:
    """Test main() orchestrates update, match-ref, workflow, and reporting on success.
    
    Verifies that:
    - The update step is run when skip_updates is False.
    - Match-ref processing is triggered when MATCH-REF column is True.
    - Snakemake workflow runs with stage='MAIN'.
    - Report is written with 'Success' status.
    - Preset warnings are shown with correct content and disabled flag.
    - System exits with code 0.
    """
    parsed_input = _build_parsed_input(skip_updates=False, disable_presets=False, match_ref=True)
    recorded: dict[str, Any] = {"updated": False, "match_ref_called": False}

    monkeypatch.setattr(vc_main, "CLIparser", lambda input_args, settings_path: parsed_input)
    monkeypatch.setattr(vc_main, "get_preset_warning_list", lambda _df: (["fallback"], ["warning"]))

    def fake_update(argv: list[str], user_config: dict[str, str]) -> None:
        """Mock update function that records invocation and validates arguments."""
        # Behavior-level check: update uses parsed configuration and is triggered when not skipped.
        assert user_config == {"profile": "x"}
        assert argv
        recorded["updated"] = True

    def fake_process_match_ref(parsed_input: Any, scheduler: Any) -> Any:
        """Mock match-ref processor that records scheduler and returns parsed input unchanged."""
        assert scheduler == "LOCAL"
        recorded["match_ref_called"] = True
        return parsed_input

    used_workflow_config = SimpleNamespace(resource_settings=object(), output_settings=object())

    def fake_run_snakemake_workflow(inputs_obj: Any, stage: str, scheduler: Any) -> tuple[bool, Any]:
        """Mock workflow executor that returns success and validates workflow parameters."""
        assert inputs_obj is parsed_input
        assert stage == "MAIN"
        assert scheduler == "LOCAL"
        recorded["run_called"] = True
        return True, used_workflow_config

    def fake_write_report(*args: Any, **kwargs: Any) -> None:
        """Mock report writer that records all report invocation arguments."""
        recorded["report"] = (args, kwargs)

    def fake_show_preset_warnings(warnings: list[str], fallbacks: list[str], disabled: bool, samples_df: pd.DataFrame) -> None:
        """Mock warning displayer that records all warning invocation arguments."""
        recorded["show"] = (warnings, fallbacks, disabled, samples_df)

    monkeypatch.setattr(vc_main, "update", fake_update)
    monkeypatch.setattr(vc_main, "process_match_ref", fake_process_match_ref)
    monkeypatch.setattr(vc_main, "run_snakemake_workflow", fake_run_snakemake_workflow)
    monkeypatch.setattr(vc_main, "WriteReport", fake_write_report)
    monkeypatch.setattr(vc_main, "show_preset_warnings", fake_show_preset_warnings)
    monkeypatch.setattr(vc_main.sys, "argv", ["viroconstrictor", "--samples", "samples.tsv"])

    with pytest.raises(SystemExit) as exc_info:
        vc_main.main(args=["--samples", "samples.tsv"], settings="/tmp/profile.ini")

    assert exc_info.value.code == 0
    assert recorded["updated"] is True
    assert recorded["match_ref_called"] is True
    assert recorded["run_called"] is True
    assert recorded["report"][0][-1] == "Success"
    assert recorded["show"][0] == ["warning"]
    assert recorded["show"][1] == ["fallback"]
    assert recorded["show"][2] is False


def test_main_failure_path_skips_update_match_ref_and_exits_one(monkeypatch: pytest.MonkeyPatch) -> None:
    """Test main() skips update/match-ref and exits with code 1 on workflow failure.
    
    Verifies that:
    - Update is skipped when skip_updates is True.
    - Match-ref is not called when not required.
    - Snakemake workflow failure results in exit code 1.
    - Report is written with 'Failed' status.
    - Preset warnings are shown with disabled flag set to True.
    """
    parsed_input = _build_parsed_input(skip_updates=True, disable_presets=True, match_ref=False)
    recorded: dict[str, Any] = {"updated": False, "match_ref_called": False}

    monkeypatch.setattr(vc_main, "CLIparser", lambda input_args, settings_path: parsed_input)
    monkeypatch.setattr(vc_main, "get_preset_warning_list", lambda _df: (["fallback"], ["warning"]))
    monkeypatch.setattr(vc_main, "update", lambda *_args, **_kwargs: recorded.__setitem__("updated", True))
    monkeypatch.setattr(vc_main, "process_match_ref", lambda *_args, **_kwargs: recorded.__setitem__("match_ref_called", True))

    used_workflow_config = SimpleNamespace(resource_settings=object(), output_settings=object())

    def fake_run_snakemake_workflow(inputs_obj: Any, stage: str, scheduler: Any) -> tuple[bool, Any]:
        """Mock workflow executor that returns failure status and validates parameters."""
        assert inputs_obj is parsed_input
        assert stage == "MAIN"
        assert scheduler == "LOCAL"
        return False, used_workflow_config

    monkeypatch.setattr(vc_main, "run_snakemake_workflow", fake_run_snakemake_workflow)

    report_calls: list[tuple[Any, ...]] = []
    monkeypatch.setattr(vc_main, "WriteReport", lambda *args, **kwargs: report_calls.append(args))

    show_calls: list[tuple[Any, Any, Any, Any]] = []
    monkeypatch.setattr(
        vc_main,
        "show_preset_warnings",
        lambda warnings, fallbacks, disabled, samples_df: show_calls.append((warnings, fallbacks, disabled, samples_df)),
    )

    with pytest.raises(SystemExit) as exc_info:
        vc_main.main()

    assert exc_info.value.code == 1
    assert recorded["updated"] is False
    assert recorded["match_ref_called"] is False
    assert report_calls[0][-1] == "Failed"
    assert show_calls[0][0] == ["warning"]
    assert show_calls[0][1] == ["fallback"]
    assert show_calls[0][2] is True
