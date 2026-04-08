"""Unit tests for :mod:`ViroConstrictor.workflow.helpers.generic_workflow_methods`.

These tests validate intended workflow-facing behavior used by the main and
match-ref Snakefiles.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import ViroConstrictor.workflow.helpers.generic_workflow_methods as methods


def _write_fasta(path: Path, records: list[tuple[str, str, str]]) -> None:
    """Write a minimal FASTA file from (id, description_suffix, sequence)."""
    lines: list[str] = []
    for rec_id, suffix, seq in records:
        if suffix:
            lines.append(f">{rec_id} {suffix}")
        else:
            lines.append(f">{rec_id}")
        lines.append(seq)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def test_get_reference_header_returns_all_record_ids(tmp_path: Path) -> None:
    """Reference headers should preserve FASTA record order and IDs."""
    fasta = tmp_path / "refs.fasta"
    _write_fasta(
        fasta,
        [
            ("RefA", "segmentA|something", "ACGT"),
            ("RefB", "segmentB|something", "TGCA"),
        ],
    )

    assert methods.get_reference_header(fasta) == ["RefA", "RefB"]


def test_read_fasta_returns_seqrecord_objects(tmp_path: Path) -> None:
    """read_fasta should parse valid FASTA input into SeqRecord objects."""
    fasta = tmp_path / "one.fasta"
    _write_fasta(fasta, [("Ref1", "segmentX|meta", "ACGTTT")])

    records = methods.read_fasta(str(fasta))

    assert len(records) == 1
    assert records[0].id == "Ref1"
    assert str(records[0].seq) == "ACGTTT"


def test_get_aminoacid_features_populates_matching_refid(monkeypatch: pytest.MonkeyPatch) -> None:
    """AA features should be set from AminoExtract output for each matching RefID."""

    monkeypatch.setattr(methods, "get_preset_parameter", lambda *_args, **_kwargs: ["CDS"])
    monkeypatch.setattr(
        methods.AminoExtract,
        "get_feature_name_attribute",
        lambda **_kwargs: {"ref-a": ["HA", "NA"], "ref-b": ["M1"]},
    )

    sample_df = pd.DataFrame(
        [
            {
                "FEATURES": "features.gff",
                "REFERENCE": "reference.fasta",
                "PRESET": "DEFAULT",
                "RefID": "ref-a",
            },
            {
                "FEATURES": "features.gff",
                "REFERENCE": "reference.fasta",
                "PRESET": "DEFAULT",
                "RefID": "ref-b",
            },
            {
                "FEATURES": "NONE",
                "REFERENCE": "reference.fasta",
                "PRESET": "DEFAULT",
                "RefID": "ref-c",
            },
        ]
    )

    out = methods.get_aminoacid_features(sample_df)

    assert out.loc[0, "AA_FEAT_NAMES"] == ("HA", "NA")
    assert out.loc[1, "AA_FEAT_NAMES"] == ("M1",)
    assert pd.isna(out.loc[2, "AA_FEAT_NAMES"])


def test_get_aminoacid_features_sets_nan_when_no_features_found(monkeypatch: pytest.MonkeyPatch) -> None:
    """If no features are returned from AminoExtract, output should be NaN."""

    monkeypatch.setattr(methods, "get_preset_parameter", lambda *_args, **_kwargs: ["CDS"])
    monkeypatch.setattr(methods.AminoExtract, "get_feature_name_attribute", lambda **_kwargs: {})

    sample_df = pd.DataFrame(
        [
            {
                "FEATURES": "features.gff",
                "REFERENCE": "reference.fasta",
                "PRESET": "DEFAULT",
                "RefID": "ref-a",
            }
        ]
    )

    out = methods.get_aminoacid_features(sample_df)

    assert pd.isna(out.loc[0, "AA_FEAT_NAMES"])


@pytest.mark.xfail(reason="Intended behavior: rows without a matching RefID should be set to NaN rather than leaving AA_FEAT_NAMES undefined.")
def test_get_aminoacid_features_unmatched_refid_defaults_to_nan(monkeypatch: pytest.MonkeyPatch) -> None:
    """Rows with no RefID match in AminoExtract output should still get NaN."""

    monkeypatch.setattr(methods, "get_preset_parameter", lambda *_args, **_kwargs: ["CDS"])
    monkeypatch.setattr(methods.AminoExtract, "get_feature_name_attribute", lambda **_kwargs: {"different-ref": ["HA"]})

    sample_df = pd.DataFrame(
        [
            {
                "FEATURES": "features.gff",
                "REFERENCE": "reference.fasta",
                "PRESET": "DEFAULT",
                "RefID": "ref-a",
            }
        ]
    )

    out = methods.get_aminoacid_features(sample_df)

    assert "AA_FEAT_NAMES" in out.columns
    assert pd.isna(out.loc[0, "AA_FEAT_NAMES"])


def test_list_aminoacid_result_outputs_builds_unique_expected_paths(monkeypatch: pytest.MonkeyPatch) -> None:
    """Generated AA result targets should be de-duplicated and correctly formatted."""

    monkeypatch.setattr(methods, "res", "results/")
    monkeypatch.setattr(methods, "amino", "aminoacids/")

    samples_df = pd.DataFrame(
        [
            {"Virus": "FluA", "RefID": "Ref1", "AA_FEAT_NAMES": ("HA", "NA")},
            {"Virus": "FluA", "RefID": "Ref1", "AA_FEAT_NAMES": ("HA",)},
            {"Virus": "FluB", "RefID": "Ref2", "AA_FEAT_NAMES": np.nan},
        ]
    )

    outputs = methods.list_aminoacid_result_outputs(samples_df)

    assert set(outputs) == {
        "results/Virus~FluA/RefID~Ref1/aminoacids/HA.faa",
        "results/Virus~FluA/RefID~Ref1/aminoacids/NA.faa",
    }


@pytest.mark.xfail(reason="Intended behavior: rows with missing/non-iterable AA feature values should be ignored without raising.")
def test_list_aminoacid_result_outputs_ignores_none_feature_values(monkeypatch: pytest.MonkeyPatch) -> None:
    """Non-list-like missing values (e.g., None) should be treated as no features."""

    monkeypatch.setattr(methods, "res", "results/")
    monkeypatch.setattr(methods, "amino", "aminoacids/")

    samples_df = pd.DataFrame([{"Virus": "FluA", "RefID": "Ref1", "AA_FEAT_NAMES": None}])

    assert methods.list_aminoacid_result_outputs(samples_df) == []


def test_segmented_ref_groups_assigns_segments_and_drops_single_segment(tmp_path: Path) -> None:
    """Segmented references need at least two unique groups; non-segmented get None."""
    multi_segment = tmp_path / "multi_segment.fasta"
    _write_fasta(
        multi_segment,
        [
            ("id1", "segA|meta", "AAAA"),
            ("id2", "segB|meta", "CCCC"),
        ],
    )

    single_segment = tmp_path / "single_segment.fasta"
    _write_fasta(
        single_segment,
        [
            ("id3", "segA|meta", "GGGG"),
            ("id4", "segA|meta", "TTTT"),
        ],
    )

    df = pd.DataFrame(
        [
            {"sample": "s1", "SEGMENTED": True, "REFERENCE": str(multi_segment)},
            {"sample": "s2", "SEGMENTED": True, "REFERENCE": str(single_segment)},
            {"sample": "s3", "SEGMENTED": False, "REFERENCE": str(single_segment)},
        ]
    )

    out = methods.segmented_ref_groups(df.copy())

    assert set(out["sample"]) == {"s1", "s3"}
    row_s1 = out[out["sample"] == "s1"].iloc[0]
    row_s3 = out[out["sample"] == "s3"].iloc[0]
    assert row_s1["segment"] == {"segA", "segB"}
    assert row_s3["segment"] == {"None"}


@pytest.mark.xfail(reason="Intended behavior: malformed FASTA descriptions should be handled gracefully, e.g. by dropping invalid rows, not by raising IndexError.")
def test_segmented_ref_groups_handles_malformed_description_without_crash(tmp_path: Path) -> None:
    """Malformed FASTA descriptions should not crash segment extraction."""
    malformed = tmp_path / "malformed.fasta"
    _write_fasta(
        malformed,
        [
            ("id1", "", "AAAA"),
            ("id2", "", "CCCC"),
        ],
    )

    df = pd.DataFrame([{"sample": "s1", "SEGMENTED": True, "REFERENCE": str(malformed)}])

    out = methods.segmented_ref_groups(df.copy())

    assert out.empty


def test_get_features_all_samples_returns_unique_features_only() -> None:
    """All-sample aggregation should deduplicate and skip missing values."""
    samples_df = pd.DataFrame(
        [
            {"AA_FEAT_NAMES": ("HA", "NA")},
            {"AA_FEAT_NAMES": ["NA", "M1"]},
            {"AA_FEAT_NAMES": np.nan},
        ]
    )

    features = methods.get_features_all_samples(samples_df)

    assert set(features) == {"HA", "NA", "M1"}


def test_get_features_per_sample_returns_unique_features_for_requested_sample() -> None:
    """Per-sample aggregation should combine all sample rows and deduplicate."""
    samples_df = pd.DataFrame(
        [
            {"sample": "S1", "AA_FEAT_NAMES": ("HA", "NA")},
            {"sample": "S1", "AA_FEAT_NAMES": ("NA", "M1")},
            {"sample": "S2", "AA_FEAT_NAMES": ("PB1",)},
        ]
    )

    assert set(methods.get_features_per_sample("S1", samples_df)) == {"HA", "NA", "M1"}


def test_get_features_per_sample_returns_empty_for_unknown_sample() -> None:
    """Unknown sample IDs should return an empty feature list."""
    samples_df = pd.DataFrame([{"sample": "S1", "AA_FEAT_NAMES": ("HA",)}])

    assert methods.get_features_per_sample("MISSING", samples_df) == []


def test_get_features_per_virus_returns_unique_features_for_requested_virus() -> None:
    """Per-virus aggregation should combine all virus rows and deduplicate."""
    samples_df = pd.DataFrame(
        [
            {"Virus": "FluA", "AA_FEAT_NAMES": ("HA", "NA")},
            {"Virus": "FluA", "AA_FEAT_NAMES": ("NA", "M1")},
            {"Virus": "FluB", "AA_FEAT_NAMES": ("PB1",)},
        ]
    )

    assert set(methods.get_features_per_virus("FluA", samples_df)) == {"HA", "NA", "M1"}


def test_get_features_per_virus_returns_empty_for_unknown_virus() -> None:
    """Unknown virus names should return an empty feature list."""
    samples_df = pd.DataFrame([{"Virus": "FluA", "AA_FEAT_NAMES": ("HA",)}])

    assert methods.get_features_per_virus("UnknownVirus", samples_df) == []


class _FakeFrame:
    """Minimal frame-like object for testing call-stack based helpers."""

    def __init__(self, f_back=None, f_locals=None, f_lineno=None) -> None:
        self.f_back = f_back
        self.f_locals = f_locals if f_locals is not None else {}
        self.f_lineno = f_lineno


def test_get_rule_name_returns_closest_rule_by_lineno(monkeypatch: pytest.MonkeyPatch) -> None:
    """Rule lookup should resolve to the closest decorator at or before call line."""
    code = "\n".join(
        [
            "@workflow.rule(name=\"prepare_refs\", lineno=10)",
            "@workflow.rule(name=\"align_reads\", lineno=40)",
            "@workflow.rule(name=\"make_consensus\", lineno=90)",
        ]
    )

    rule_frame = _FakeFrame(f_locals={"code": code})
    caller = _FakeFrame(f_back=rule_frame, f_lineno=57)
    frame = _FakeFrame(f_back=caller)

    monkeypatch.setattr(methods.inspect, "currentframe", lambda: frame)

    assert methods.get_rule_name() == "align_reads"


def test_get_rule_name_returns_fallback_when_no_rule_matches(monkeypatch: pytest.MonkeyPatch) -> None:
    """If all discovered rules are after the call line, a safe fallback is expected."""
    code = "\n".join(
        [
            "@workflow.rule(name=\"late_rule\", lineno=100)",
            "@workflow.rule(name=\"later_rule\", lineno=150)",
        ]
    )

    rule_frame = _FakeFrame(f_locals={"code": code})
    caller = _FakeFrame(f_back=rule_frame, f_lineno=25)
    frame = _FakeFrame(f_back=caller)

    monkeypatch.setattr(methods.inspect, "currentframe", lambda: frame)

    assert methods.get_rule_name() == "unknown_rule"


def test_get_rule_name_returns_fallback_with_missing_stack_context(monkeypatch: pytest.MonkeyPatch) -> None:
    """Missing frame/call-site metadata should not crash rule-name resolution."""
    monkeypatch.setattr(methods.inspect, "currentframe", lambda: None)

    assert methods.get_rule_name() == "unknown_rule"


def test_get_rule_name_returns_fallback_when_code_not_available(monkeypatch: pytest.MonkeyPatch) -> None:
    """When calling context does not provide workflow source code, fallback is used."""
    rule_frame = _FakeFrame(f_locals={"not_code": "..."})
    caller = _FakeFrame(f_back=rule_frame, f_lineno=42)
    frame = _FakeFrame(f_back=caller)

    monkeypatch.setattr(methods.inspect, "currentframe", lambda: frame)

    assert methods.get_rule_name() == "unknown_rule"
