import sys
from argparse import ArgumentParser
from pathlib import Path

import pandas as pd
import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(PROJECT_ROOT.joinpath("ViroConstrictor/workflow")))
from ViroConstrictor.workflow.main.scripts.group_aminoacids import GroupAminoAcids  # isort:skip

def _create_space_pickle(path: Path) -> None:
    pd.DataFrame(
        {
            "sample": ["sample1", "sample2"],
            "Virus": ["VirusA", "VirusA"],
            "RefID": ["RefA", "RefA"],
            "AA_FEAT_NAMES": [["S", "N"], ["S", "N"]],
        }
    ).to_pickle(path)


def test_add_arguments_parses_space_path() -> None:
    parser = ArgumentParser()
    GroupAminoAcids.add_arguments(parser)
    args = parser.parse_args(["--input", "in.faa", "--output", "out.faa", "--space", "space.pkl"])
    assert str(args.space).endswith("space.pkl")


def test_group_amino_acids_writes_feature_grouped_outputs(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.chdir(tmp_path)
    input_fasta = tmp_path / "input.faa"
    space_pkl = tmp_path / "space.pkl"
    _create_space_pickle(space_pkl)

    input_fasta.write_text(
        ">sample1.S\nAAAA\n"
        ">sample1.N\nTTTT\n"
        ">sample2.S\nCCCC\n",
        encoding="utf-8",
    )

    out_s = "results/Virus~VirusA/RefID~RefA/aminoacids/S.faa"
    out_n = "results/Virus~VirusA/RefID~RefA/aminoacids/N.faa"
    output_targets = f"{out_s} {out_n}"

    GroupAminoAcids(input=str(input_fasta), output=output_targets, space=str(space_pkl)).run()

    s_content = (tmp_path / out_s).read_text(encoding="utf-8")
    n_content = (tmp_path / out_n).read_text(encoding="utf-8")
    assert ">sample1.S" in s_content
    assert ">sample2.S" in s_content
    assert ">sample1.N" in n_content


def test_group_amino_acids_writes_only_requested_output_targets(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.chdir(tmp_path)
    input_fasta = tmp_path / "input.faa"
    space_pkl = tmp_path / "space.pkl"
    _create_space_pickle(space_pkl)

    input_fasta.write_text(
        ">sample1.S\nAAAA\n"
        ">sample1.N\nTTTT\n",
        encoding="utf-8",
    )

    out_s = "results/Virus~VirusA/RefID~RefA/aminoacids/S.faa"
    out_n = "results/Virus~VirusA/RefID~RefA/aminoacids/N.faa"

    GroupAminoAcids(input=str(input_fasta), output=out_s, space=str(space_pkl)).run()

    assert (tmp_path / out_s).exists()
    assert not (tmp_path / out_n).exists()


def test_group_amino_acids_requires_expected_space_columns(tmp_path: Path) -> None:
    input_fasta = tmp_path / "input.faa"
    space_pkl = tmp_path / "space.pkl"
    input_fasta.write_text(">sample1.S\nAAAA\n", encoding="utf-8")

    pd.DataFrame({"sample": ["sample1"]}).to_pickle(space_pkl)

    with pytest.raises(KeyError):
        GroupAminoAcids(input=str(input_fasta), output="results/Virus~V/RefID~R/aminoacids/S.faa", space=str(space_pkl)).run()


@pytest.mark.xfail(reason="Empty parsed input currently crashes on missing Virus column; intended behavior is graceful no-op output", strict=False)
def test_group_amino_acids_handles_empty_input_without_crashing(tmp_path: Path) -> None:
    input_fasta = tmp_path / "empty.faa"
    space_pkl = tmp_path / "space.pkl"
    input_fasta.write_text("", encoding="utf-8")
    _create_space_pickle(space_pkl)

    GroupAminoAcids(
        input=str(input_fasta),
        output="results/Virus~VirusA/RefID~RefA/aminoacids/S.faa",
        space=str(space_pkl),
    ).run()
