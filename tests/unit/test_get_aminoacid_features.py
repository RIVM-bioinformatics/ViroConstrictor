"""Unit tests for get_aminoacid_features helper."""

import numpy as np
import pandas as pd

from ViroConstrictor.workflow.helpers import generic_workflow_methods as gwm


def test_get_aminoacid_features_sets_nan_when_features_none() -> None:
    """Rows with FEATURES='NONE' should have NaN amino acid feature names."""
    df = pd.DataFrame(
        [
            {
                "FEATURES": "NONE",
                "REFERENCE": "ref.fa",
                "PRESET": "SARSCOV2",
                "RefID": "REF_A",
            }
        ]
    )

    out = gwm.get_aminoacid_features(df)

    assert "AA_FEAT_NAMES" in out.columns
    assert np.isnan(out.loc[0, "AA_FEAT_NAMES"])


def test_get_aminoacid_features_assigns_tuple_for_matching_refid(monkeypatch) -> None:
    """Feature names should be assigned as tuples when RefID is found in AminoExtract output."""

    def _fake_feature_types(_preset_name: str) -> list[str]:
        return ["gene", "CDS"]

    def _fake_get_feature_name_attribute(input_gff: str, input_seq: str, feature_types: list[str]):
        assert input_gff == "features.gff"
        assert input_seq == "reference.fa"
        assert feature_types == ["gene", "CDS"]
        return {
            "REF_A": ["S", "N"],
            "REF_OTHER": ["M"],
        }

    monkeypatch.setattr(gwm, "get_feature_types_from_preset", _fake_feature_types)
    monkeypatch.setattr(gwm.AminoExtract, "get_feature_name_attribute", _fake_get_feature_name_attribute)

    df = pd.DataFrame(
        [
            {
                "FEATURES": "features.gff",
                "REFERENCE": "reference.fa",
                "PRESET": "SARSCOV2",
                "RefID": "REF_A",
            }
        ]
    )

    out = gwm.get_aminoacid_features(df)

    assert out.loc[0, "AA_FEAT_NAMES"] == ("S", "N")


def test_get_feature_types_from_preset_parses_singular_option(monkeypatch) -> None:
    """--feature-type should be parsed into a list of feature names."""

    monkeypatch.setattr(
        gwm,
        "get_flags",
        lambda _stage, _preset, _task, _platform: "--keep-gaps --feature-type gene CDS",
    )

    assert gwm.get_feature_types_from_preset("SARSCOV2") == ["gene", "CDS"]


def test_get_feature_types_from_preset_parses_plural_option(monkeypatch) -> None:
    """--feature-types should also be parsed for backward compatibility."""

    monkeypatch.setattr(
        gwm,
        "get_flags",
        lambda _stage, _preset, _task, _platform: "--keep-gaps --feature-types gene",
    )

    assert gwm.get_feature_types_from_preset("ORTHOHEPADNAVIRUS") == ["gene"]
