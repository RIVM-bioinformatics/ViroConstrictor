"""Unit tests for :mod:`ViroConstrictor.workflow.helpers.presets`.

Tests cover fuzzy matching, explicit and implicit stage resolution, fallback to
DEFAULT preset values, and failure handling.
"""

from __future__ import annotations

import pytest

import ViroConstrictor.workflow.helpers.presets as presets


def test_get_key_from_value_returns_matching_group_name_or_none() -> None:
    """Find the preset group that contains a provided alias."""
    aliases = {
        "GROUP_A": ["A1", "A2"],
        "GROUP_B": ["B1"],
    }

    assert presets.get_key_from_value(aliases, "A2") == "GROUP_A"
    assert presets.get_key_from_value(aliases, "missing") is None


def test_match_preset_name_returns_default_when_presets_disabled() -> None:
    """When presets are disabled, matching should not be attempted."""
    matched_name, score = presets.match_preset_name("SARS-CoV-2", use_presets=False)
    assert matched_name == "DEFAULT"
    assert score == 0.0


def test_match_preset_name_normalizes_input_and_matches_alias() -> None:
    """Special characters should be stripped before matching aliases."""
    matched_name, score = presets.match_preset_name("sars-cov-2!!!", use_presets=True)
    assert matched_name == "SARSCOV2"
    assert score >= 0.40


def test_match_preset_name_respects_explicit_default_query() -> None:
    """DEFAULT should map directly to DEFAULT with full confidence."""
    matched_name, score = presets.match_preset_name("DEFAULT", use_presets=True)
    assert matched_name == "DEFAULT"
    assert score == 1.0


def test_match_preset_name_falls_back_to_default_for_low_similarity() -> None:
    """Unknown targets with poor similarity should not map to a specific preset."""
    matched_name, score = presets.match_preset_name("zzzzzzzzzzzz", use_presets=True)
    assert matched_name == "DEFAULT"
    assert score == 0.0


@pytest.mark.xfail(reason="Intended behavior: empty alias tables should safely return DEFAULT instead of raising IndexError.")
def test_match_preset_name_handles_empty_alias_table_gracefully(monkeypatch: pytest.MonkeyPatch) -> None:
    """Matching should remain safe even when alias data cannot be loaded."""
    monkeypatch.setattr(presets, "aliases", {})
    matched_name, score = presets.match_preset_name("anything", use_presets=True)
    assert matched_name == "DEFAULT"
    assert score == 0.0


def test_collapse_preset_group_collects_selected_stages_with_prefix(monkeypatch: pytest.MonkeyPatch) -> None:
    """Selected stage keys should be flattened with the stage identifier prefix."""
    fake_presets = {
        "CUSTOM": {
            "STAGE_MAIN": {"A": "x", "B": "y"},
            "STAGE_MATCHREF": {"C": "z"},
            "IGNORED": {"D": "q"},
        }
    }
    monkeypatch.setattr(presets, "presets", fake_presets)

    collapsed = presets.collapse_preset_group("CUSTOM", ["STAGE_MAIN", "STAGE_MATCHREF"], "STAGE_MAIN")

    assert collapsed == {
        "STAGE_MAIN_A": "x",
        "STAGE_MAIN_B": "y",
        "STAGE_MAIN_C": "z",
    }


def test_get_preset_parameter_returns_preset_override_when_present() -> None:
    """Preset-specific values should override DEFAULT for matching stage."""
    value = presets.get_preset_parameter(
        preset_name="INFLUENZA",
        parameter_name="Minimap2_Settings",
        stage_identifier="STAGE_MAIN",
    )

    assert value == "--frag=no"


def test_get_preset_parameter_falls_back_to_default_when_missing_from_preset() -> None:
    """Missing stage parameter in a preset should inherit from DEFAULT."""
    value = presets.get_preset_parameter(
        preset_name="INFLUENZA",
        parameter_name="AmpliGone_AlignmentPreset_illumina",
        stage_identifier="STAGE_MAIN",
    )

    assert value == "-ap sr"


def test_get_preset_parameter_supports_matchref_stage_values() -> None:
    """Match-ref stage requests should return match-ref-specific settings."""
    value = presets.get_preset_parameter(
        preset_name="INFLUENZA",
        parameter_name="Samtools_Filters_nanopore",
        stage_identifier="STAGE_MATCHREF",
    )

    assert value == "-q 5"


def test_get_preset_parameter_uses_callers_vc_stage_when_not_passed() -> None:
    """If stage_identifier is omitted, caller global VC_STAGE should be used."""
    # Required by presets.get_preset_parameter inspect-based stage lookup.
    global VC_STAGE
    VC_STAGE = "STAGE_MAIN"

    def _call_without_explicit_stage() -> str:
        """Call get_preset_parameter without explicit stage_identifier argument.

        Returns
        -------
        str
            Preset parameter value retrieved using caller's VC_STAGE global.
        """
        return presets.get_preset_parameter(
            preset_name="DEFAULT",
            parameter_name="Fastp_MinReadLength_Settings_illumina",
        )

    assert _call_without_explicit_stage() == "--length_required 100"


def test_get_preset_parameter_raises_for_unknown_preset() -> None:
    """Unknown preset names should raise a clear key lookup failure."""
    with pytest.raises(KeyError):
        presets.get_preset_parameter(
            preset_name="DOES_NOT_EXIST",
            parameter_name="Minimap2_Settings",
            stage_identifier="STAGE_MAIN",
        )


def test_get_preset_parameter_raises_when_parameter_missing_in_default(monkeypatch: pytest.MonkeyPatch) -> None:
    """If no preset and no DEFAULT value exist, the function should raise."""
    fake_presets = {
        "DEFAULT": {
            "STAGE_MAIN": {},
            "STAGE_MATCHREF": {},
            "STAGE_GLOBAL": {},
        },
        "CUSTOM": {
            "STAGE_MAIN": {},
            "STAGE_MATCHREF": {},
            "STAGE_GLOBAL": {},
        },
    }
    monkeypatch.setattr(presets, "presets", fake_presets)

    with pytest.raises(KeyError):
        presets.get_preset_parameter(
            preset_name="CUSTOM",
            parameter_name="Missing_Parameter",
            stage_identifier="STAGE_MAIN",
        )
