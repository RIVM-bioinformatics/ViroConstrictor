"""
Unit tests for CombineTabular script.

Tests concatenation of tabular result files (coverage, mutations, amplicon
coverage) with per-sample metadata enrichment and file-type-specific formatting.
"""

import sys
from argparse import ArgumentParser
from pathlib import Path

import pandas as pd
import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(PROJECT_ROOT.joinpath("ViroConstrictor/workflow")))
from ViroConstrictor.workflow.main.scripts.combine_tabular import CombineTabular  # isort:skip


def _combiner(
    *,
    output: Path,
    input_files: list[Path | str],
    virus_list: list[str],
    refid_list: list[str],
    file_type: str,
    separator: str = "\t",
) -> CombineTabular:
    """
    Factory function to create CombineTabular instance with common test parameters.

    Parameters
    ----------
    output : Path
        Output file path.
    input_files : list[Path | str]
        Input tabular files to combine.
    virus_list : list[str]
        Virus names corresponding to each input file.
    refid_list : list[str]
        Reference IDs corresponding to each input file.
    file_type : str
        File type identifier ("coverage", "mutations", "amplicon_coverage").
    separator : str, optional
        Field separator (tab or comma). Default: "\t".

    Returns
    -------
    CombineTabular
        Configured CombineTabular instance ready for testing.
    """
    return CombineTabular(
        input="",
        output=output,
        input_files=input_files,
        virus_list=virus_list,
        refid_list=refid_list,
        file_type=file_type,  # type: ignore[arg-type]
        separator=separator,
    )


def test_add_arguments_registers_expected_flags() -> None:
    """
    Verify that CombineTabular CLI parser registers all expected flags.

    Tests that the argument parser accepts input_files, output, virus_list,
    refid_list, file_type, and separator flags and correctly stores their values.
    """
    parser = ArgumentParser()
    CombineTabular.add_arguments(parser)
    args = parser.parse_args(
        [
            "--input",
            "dummy",
            "--output",
            "out.tsv",
            "--input_files",
            "a.tsv",
            "b.tsv",
            "--virus_list",
            "VirusA",
            "VirusB",
            "--refid_list",
            "RefA",
            "RefB",
            "--file_type",
            "coverage",
        ]
    )
    assert args.input_files == ["a.tsv", "b.tsv"]
    assert args.virus_list == ["VirusA", "VirusB"]
    assert args.refid_list == ["RefA", "RefB"]
    assert args.file_type == "coverage"
    assert args.separator == "\t"


def test_iter_nonempty_files_excludes_missing_and_empty(tmp_path: Path) -> None:
    """
    Verify that empty and missing files are excluded from iteration.

    Tests that _iter_nonempty_files() correctly filters out files that do not
    exist or are empty, returning only non-empty existing files.
    """
    non_empty = tmp_path / "has_content.tsv"
    empty = tmp_path / "empty.tsv"
    missing = tmp_path / "missing.tsv"
    non_empty.write_text("x\n")
    empty.write_text("")

    combiner = _combiner(
        output=tmp_path / "out.tsv",
        input_files=[missing, empty, non_empty],
        virus_list=["v1", "v2", "v3"],
        refid_list=["r1", "r2", "r3"],
        file_type="coverage",
    )

    assert list(combiner._iter_nonempty_files()) == [non_empty]


def test_read_tabular_file_returns_none_for_empty_file(tmp_path: Path) -> None:
    """
    Verify that reading an empty tabular file returns None.

    Tests that _read_tabular_file() handles empty files gracefully by returning None
    instead of an empty DataFrame.
    """
    infile = tmp_path / "empty.tsv"
    infile.write_text("")
    combiner = _combiner(
        output=tmp_path / "out.tsv",
        input_files=[infile],
        virus_list=["v"],
        refid_list=["r"],
        file_type="mutations",
    )

    assert combiner._read_tabular_file(infile) is None


def test_drop_header_row_if_present_only_removes_true_header() -> None:
    """
    Verify that header row is only removed when it matches expected header.

    Tests that _drop_header_row_if_present() removes a matching header row but
    leaves non-matching rows untouched.
    """
    expected = ["A", "B"]
    df_with_header = pd.DataFrame([expected, ["1", "2"]])
    df_without_header = pd.DataFrame([["x", "y"], ["1", "2"]])

    stripped = CombineTabular._drop_header_row_if_present(df_with_header, expected)
    untouched = CombineTabular._drop_header_row_if_present(df_without_header, expected)

    assert stripped.shape == (1, 2)
    assert stripped.iloc[0].tolist() == ["1", "2"]
    assert untouched.equals(df_without_header)


def test_combine_coverage_merges_valid_files_and_adds_metadata(tmp_path: Path) -> None:
    """
    Verify coverage files are merged with virus and reference ID metadata columns.

    Tests that the _combine_coverage() method correctly merges coverage data from
    multiple files while adding Virus and Reference_ID columns in the correct order.
    """
    f1 = tmp_path / "coverage1.tsv"
    f2 = tmp_path / "coverage2.tsv"
    out = tmp_path / "combined_coverage.tsv"

    f1.write_text(
        "Sample_name\tWidth_at_mincov_1\tWidth_at_mincov_5\tWidth_at_mincov_10\tWidth_at_mincov_50\tWidth_at_mincov_100\n"
        "sample1\t10\t20\t30\t40\t50\n"
    )
    f2.write_text("sample2\t11\t21\t31\t41\t51\n")

    combiner = _combiner(
        output=out,
        input_files=[f1, f2],
        virus_list=["VirusA", "VirusB"],
        refid_list=["NA", "RefB"],
        file_type="coverage",
    )

    combiner.run()
    df = pd.read_csv(out, sep="\t", keep_default_na=False)

    assert list(df.columns) == [
        "Sample_name",
        "Virus",
        "Reference_ID",
        "Width_at_mincov_1",
        "Width_at_mincov_5",
        "Width_at_mincov_10",
        "Width_at_mincov_50",
        "Width_at_mincov_100",
    ]
    assert set(df["Sample_name"]) == {"sample1", "sample2"}
    assert set(df["Virus"]) == {"VirusA", "VirusB"}
    assert df.loc[df["Sample_name"] == "sample1", "Reference_ID"].iloc[0] == "NA"


def test_combine_coverage_empty_inputs_writes_header_only(tmp_path: Path) -> None:
    """
    Verify that empty coverage inputs produce header-only output.

    Tests that when all input coverage files are empty, the output contains
    only the header row with no data rows.
    """
    empty = tmp_path / "empty.tsv"
    out = tmp_path / "combined_coverage.tsv"
    empty.write_text("")
    combiner = _combiner(
        output=out,
        input_files=[empty],
        virus_list=["VirusA"],
        refid_list=["RefA"],
        file_type="coverage",
    )

    combiner.run()
    lines = out.read_text().splitlines()
    assert len(lines) == 1
    assert lines[0].startswith("Sample_name\tVirus\tReference_ID")


def test_combine_mutations_uses_mapped_reference_and_column_order(tmp_path: Path) -> None:
    """
    Verify mutation files are combined with mapped reference IDs and correct column order.

    Tests that the _combine_mutations() method replaces the Reference_ID column
    from the file with the mapped value and maintains the expected column order.
    """
    f1 = tmp_path / "mut1.tsv"
    f2 = tmp_path / "mut2.tsv"
    out = tmp_path / "combined_mutations.tsv"

    f1.write_text("Sample\tReference_ID\tPosition\tReference_Base\tVariant_Base\tDepth\nsample1\tWRONG_FROM_FILE\t100\tA\tT\t50\n")
    f2.write_text("sample2\tIGNORED\t200\tC\tG\t60\n")

    combiner = _combiner(
        output=out,
        input_files=[f1, f2],
        virus_list=["VirusA", "VirusB"],
        refid_list=["NA", "RefB"],
        file_type="mutations",
    )

    combiner.run()
    df = pd.read_csv(out, sep="\t", keep_default_na=False)

    assert list(df.columns) == [
        "Sample",
        "Virus",
        "Reference_ID",
        "Position",
        "Reference_Base",
        "Variant_Base",
        "Depth",
    ]
    assert set(df["Sample"]) == {"sample1", "sample2"}
    assert df.loc[df["Sample"] == "sample1", "Reference_ID"].iloc[0] == "NA"
    assert df.loc[df["Sample"] == "sample2", "Reference_ID"].iloc[0] == "RefB"


def test_combine_mutations_empty_inputs_writes_header_only(tmp_path: Path) -> None:
    """
    Verify that empty mutation inputs produce header-only output.

    Tests that when all input mutation files are empty, the output contains
    only the header row with correct column names.
    """
    empty = tmp_path / "empty.tsv"
    out = tmp_path / "combined_mutations.tsv"
    empty.write_text("")
    combiner = _combiner(
        output=out,
        input_files=[empty],
        virus_list=["VirusA"],
        refid_list=["RefA"],
        file_type="mutations",
    )

    combiner.run()
    lines = out.read_text().splitlines()
    assert len(lines) == 1
    assert lines[0] == "Sample\tVirus\tReference_ID\tPosition\tReference_Base\tVariant_Base\tDepth"


def test_combine_amplicon_coverage_adds_front_metadata_columns(tmp_path: Path) -> None:
    """
    Verify that amplicon coverage files have metadata columns added to the front.

    Tests that the _combine_amplicon_coverage() method prepends Virus and Reference_ID
    columns to the front of the combined output.
    """
    f1 = tmp_path / "amp1.csv"
    f2 = tmp_path / "amp2.csv"
    out = tmp_path / "combined_amplicon.csv"

    f1.write_text("amplicon,mean_cov,median_cov\namp1,100,95\n")
    f2.write_text("amplicon,mean_cov,median_cov\namp2,150,145\n")

    combiner = _combiner(
        output=out,
        input_files=[f1, f2],
        virus_list=["VirusA", "VirusB"],
        refid_list=["NA", "RefB"],
        file_type="amplicon_coverage",
        separator=",",
    )

    combiner.run()
    df = pd.read_csv(out, keep_default_na=False)

    assert list(df.columns[:2]) == ["Virus", "Reference_ID"]
    assert set(df["amplicon"]) == {"amp1", "amp2"}
    assert df.loc[df["amplicon"] == "amp1", "Reference_ID"].iloc[0] == "NA"


def test_combine_amplicon_coverage_empty_inputs_writes_empty_csv(tmp_path: Path) -> None:
    """
    Verify that empty amplicon coverage inputs produce empty CSV with newline.

    Tests that when all input amplicon coverage files are empty, the output
    contains only a newline character.
    """
    empty = tmp_path / "empty.csv"
    out = tmp_path / "combined_amplicon.csv"
    empty.write_text("")

    combiner = _combiner(
        output=out,
        input_files=[empty],
        virus_list=["VirusA"],
        refid_list=["RefA"],
        file_type="amplicon_coverage",
        separator=",",
    )
    combiner.run()

    assert out.read_text() == "\n"


def test_run_dispatches_to_expected_handler(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """
    Verify run() method dispatches to the correct handler based on file_type.

    Tests that the run() method calls the appropriate combining method
    (_combine_coverage, _combine_mutations, or _combine_amplicon_coverage) based
    on the file_type parameter.
    """
    out = tmp_path / "out.tsv"
    infile = tmp_path / "in.tsv"
    infile.write_text("sample\t1\t2\t3\t4\t5\n")
    combiner = _combiner(
        output=out,
        input_files=[infile],
        virus_list=["VirusA"],
        refid_list=["RefA"],
        file_type="coverage",
    )

    called = {"coverage": 0}

    def fake_combine_coverage(file_map: dict[Path | str, tuple[str, str]]) -> None:
        """Mock coverage combining function for testing."""
        assert infile in file_map
        called["coverage"] += 1

    monkeypatch.setattr(combiner, "_combine_coverage", fake_combine_coverage)
    combiner.run()

    assert called["coverage"] == 1


@pytest.mark.xfail(
    strict=True,
    reason="Docstring says invalid file_type should raise ValueError, but current implementation silently no-ops.",
)
def test_invalid_file_type_raises_value_error(tmp_path: Path) -> None:
    """
    Test unsupported file type rejection (XFAIL - intended behavior).

    Intended behavior: Unsupported file_type values should raise ValueError, but
    the current implementation silently skips processing instead.
    """
    out = tmp_path / "out.tsv"
    infile = tmp_path / "in.tsv"
    infile.write_text("irrelevant\n")

    combiner = _combiner(
        output=out,
        input_files=[infile],
        virus_list=["VirusA"],
        refid_list=["RefA"],
        file_type="not_supported",
    )

    with pytest.raises(ValueError):
        combiner.run()
