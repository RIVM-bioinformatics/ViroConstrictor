"""Unit tests for main-workflow samples_df creation helper."""

from pathlib import Path

import pandas as pd

from ViroConstrictor.workflow.helpers import generic_workflow_methods as gwm


def test_create_samples_df_builds_expected_shape_and_columns(tmp_path: Path) -> None:
    """samples_df should normalize sample names, rename VIRUS, and explode reference IDs."""
    ref_path = tmp_path / "reference.fasta"
    ref_path.write_text(
        ">REF_A\nACGTACGT\n>REF_B\nTTTTGGGG\n",
        encoding="utf-8",
    )

    samples = {
        "sample_1": {
            "VIRUS": "sars",
            "REFERENCE": str(ref_path),
            "FEATURES": "NONE",
            "PRESET": "SARSCOV2",
        },
        "sample_2": {
            "VIRUS": "sars",
            "REFERENCE": str(ref_path),
            "FEATURES": "NONE",
            "PRESET": "SARSCOV2",
        },
    }

    samples_df = gwm.create_samples_df(samples)

    assert len(samples_df) == 4
    assert set(samples_df["sample"]) == {"sample_1", "sample_2"}
    assert set(samples_df["Virus"]) == {"sars"}
    assert set(samples_df["RefID"]) == {"REF_A", "REF_B"}
    assert "AA_FEAT_NAMES" in samples_df.columns
    assert samples_df["AA_FEAT_NAMES"].isna().all()
    assert all(pd.isna(value) for value in samples_df["AA_FEAT_NAMES"])
