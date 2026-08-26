"""Focused tests for presets2 helper resolution used by Snakemake rules."""

import pytest

from ViroConstrictor.workflow.helpers import presets as presets_module

# Snakemake-style global used by presets2 inspect fallback when stage_identifier is omitted.
VC_STAGE = "MAIN"


def _smk_get_flags(preset_name: str, task: str, platform: str) -> str:
    return presets_module.get_flags(
        stage_identifier="",
        preset_name=preset_name,
        task=task,
        platform=platform,
    )


@pytest.fixture
def mock_structured_presets(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(
        presets_module,
        "main_presets",
        {
            "DEFAULT": {
                "qc": {
                    "bin_by_platform": {
                        "illumina": "fastp",
                    },
                    "flags_by_platform": {
                        "illumina": "--length_required 100",
                    },
                },
                "alignment": {
                    "bin_by_platform": {
                        "nanopore": "minimap2",
                    },
                    "flags_by_platform": {
                        "nanopore": "--secondary=no",
                    },
                },
            },
            "CUSTOM": {
                "qc": {
                    "flags_by_platform": {
                        "illumina": "--length_required 150",
                    }
                }
            },
        },
    )

    monkeypatch.setattr(
        presets_module,
        "matchref_presets",
        {
            "DEFAULT": {
                "alignment": {
                    "bin_by_platform": {"nanopore": "minimap2"},
                    "flags_by_platform": {"nanopore": "--secondary=no"},
                }
            }
        },
    )


@pytest.fixture
def mock_aliases(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(
        presets_module,
        "aliases",
        {
            "SARSCOV2": ["SARS-COV-2", "SARS"],
            "DEFAULT": ["DEFAULT"],
        },
    )


def test_get_flags_prefers_preset_specific_values(mock_structured_presets: None) -> None:
    assert (
        presets_module.get_flags(
            stage_identifier="MAIN",
            preset_name="CUSTOM",
            task="qc",
            platform="illumina",
        )
        == "--length_required 150"
    )


def test_get_binary_falls_back_to_default_for_missing_preset_key(mock_structured_presets: None) -> None:
    assert (
        presets_module.get_binary(
            stage_identifier="MAIN",
            preset_name="CUSTOM",
            task="qc",
            platform="illumina",
        )
        == "fastp"
    )


def test_get_flags_can_resolve_stage_from_snakemake_style_global(mock_structured_presets: None) -> None:
    assert _smk_get_flags("CUSTOM", "qc", "illumina") == "--length_required 150"


def test_match_preset_name_returns_default_when_presets_disabled(mock_aliases: None) -> None:
    assert presets_module.match_preset_name("sars", use_presets=False) == ("DEFAULT", 0.0)


def test_match_preset_name_matches_alias_with_sufficient_similarity(mock_aliases: None) -> None:
    preset, score = presets_module.match_preset_name("sars-cov-2", use_presets=True)
    assert preset == "SARSCOV2"
    assert score >= 0.40


def test_match_preset_name_returns_default_when_similarity_too_low(mock_aliases: None) -> None:
    assert presets_module.match_preset_name("totally-unrelated", use_presets=True) == ("DEFAULT", 0.0)


def test_get_flags_raises_for_unknown_task_and_platform(mock_structured_presets: None) -> None:
    with pytest.raises(KeyError, match="Flags not found"):
        presets_module.get_flags(
            stage_identifier="MAIN",
            preset_name="CUSTOM",
            task="unknown_task",
            platform="nanopore",
        )
