import sys
import types
from argparse import ArgumentParser
from pathlib import Path

import pandas as pd
import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(PROJECT_ROOT.joinpath("ViroConstrictor/workflow")))

# Keep tests hermetic even when AminoExtract is unavailable in the environment.
if "AminoExtract" not in sys.modules:
    aminoextract_stub = types.ModuleType("AminoExtract")
    aminoextract_stub.SequenceReader = object
    sys.modules["AminoExtract"] = aminoextract_stub

from ViroConstrictor.workflow.main.scripts.extract_gff import ExtractGff  # isort:skip


class _FakeGff:
    def __init__(self, df: pd.DataFrame) -> None:
        self.df = df
        self.exported_to: Path | str | None = None

    def export_gff_to_file(self, output: Path | str) -> None:
        self.exported_to = output
        Path(output).write_text(self.df.to_csv(index=False), encoding="utf-8")


def test_add_arguments_parses_ref_id() -> None:
    parser = ArgumentParser()
    ExtractGff.add_arguments(parser)
    args = parser.parse_args(["--input", "in.gff", "--output", "out.gff", "--ref_id", "REF_1"])
    assert args.ref_id == "REF_1"


def test_extract_gff_filters_to_requested_reference(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    input_gff = tmp_path / "input.gff"
    output_gff = tmp_path / "output.gff"
    input_gff.write_text("dummy", encoding="utf-8")

    fake_gff = _FakeGff(
        pd.DataFrame(
            {
                "seqid": ["REF_A", "REF_B", "REF_A"],
                "type": ["gene", "gene", "CDS"],
                "start": [1, 5, 8],
            }
        )
    )

    class _FakeSequenceReader:
        def __init__(self, logger: object | None, verbose: bool) -> None:
            assert logger is None
            assert not verbose

        def read_gff(self, path: Path | str) -> _FakeGff:
            assert Path(path) == input_gff
            return fake_gff

    module = sys.modules[ExtractGff.__module__]
    monkeypatch.setattr(module, "SequenceReader", _FakeSequenceReader)

    extractor = ExtractGff(input=input_gff, output=output_gff, ref_id="REF_A")
    extractor.run()

    assert fake_gff.exported_to == output_gff
    assert set(fake_gff.df["seqid"]) == {"REF_A"}
    assert len(fake_gff.df) == 2


def test_extract_gff_with_no_matching_reference_writes_empty_result(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    input_gff = tmp_path / "input.gff"
    output_gff = tmp_path / "output.gff"
    input_gff.write_text("dummy", encoding="utf-8")

    fake_gff = _FakeGff(pd.DataFrame({"seqid": ["REF_A"], "type": ["gene"]}))

    class _FakeSequenceReader:
        def __init__(self, logger: object | None, verbose: bool) -> None:
            pass

        def read_gff(self, path: Path | str) -> _FakeGff:
            return fake_gff

    module = sys.modules[ExtractGff.__module__]
    monkeypatch.setattr(module, "SequenceReader", _FakeSequenceReader)

    extractor = ExtractGff(input=input_gff, output=output_gff, ref_id="DOES_NOT_EXIST")
    extractor.run()

    assert output_gff.exists()
    assert fake_gff.df.empty


def test_extract_gff_surfaces_reader_errors(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    input_gff = tmp_path / "input.gff"
    output_gff = tmp_path / "output.gff"
    input_gff.write_text("dummy", encoding="utf-8")

    class _FailingSequenceReader:
        def __init__(self, logger: object | None, verbose: bool) -> None:
            pass

        def read_gff(self, path: Path | str) -> _FakeGff:
            raise RuntimeError("malformed gff")

    module = sys.modules[ExtractGff.__module__]
    monkeypatch.setattr(module, "SequenceReader", _FailingSequenceReader)

    extractor = ExtractGff(input=input_gff, output=output_gff, ref_id="REF_A")
    with pytest.raises(RuntimeError, match="malformed gff"):
        extractor.run()


@pytest.mark.xfail(reason="Current implementation uses assert for type checks instead of explicit ValueError validation", strict=False)
def test_extract_gff_rejects_non_string_ref_id_with_clear_error(tmp_path: Path) -> None:
    input_gff = tmp_path / "input.gff"
    output_gff = tmp_path / "output.gff"
    input_gff.write_text("dummy", encoding="utf-8")

    extractor = ExtractGff(input=input_gff, output=output_gff, ref_id=123)  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="ref_id|string"):
        extractor.run()
