"""Test parser utilities and CLI parser validation paths.

These tests cover samplesheet handling, path normalization, file helpers, and
CLI parser helpers with mocked filesystem and dataframe dependencies.
"""

import tempfile
import uuid
from argparse import Namespace
from configparser import ConfigParser
from pathlib import Path

import pandas as pd
import pytest

from ViroConstrictor import scheduler as scheduler_module
from ViroConstrictor.parser import (
    CheckInputFiles,
    CLIparser,
    args_to_df,
    check_file_extension,
    check_samplesheet_columns,
    check_samplesheet_empty_rows,
    check_samplesheet_rows,
    convert_log_text,
    dir_path,
    file_exists,
    is_csv_file,
    is_excel_file,
    is_tsv_file,
    open_sample_sheet,
    required_cols,
    sampledir_to_df,
    samplesheet_enforce_absolute_paths,
)


def _build_args(tmp_path: Path, **overrides) -> Namespace:
    """Build a parser-like namespace for parser helper tests.

    Parameters
    ----------
    tmp_path : Path
        Temporary directory used to create the reference file.
    **overrides : Any
        Keyword overrides for the default CLI argument values.

    Returns
    -------
    Namespace
        Namespace populated with default parser inputs.
    """
    reference = tmp_path / "ref.fasta"
    reference.write_text(">ref\nACGT\n", encoding="utf-8")
    defaults = {
        "input": str(tmp_path),
        "platform": "nanopore",
        "reference": str(reference),
        "primers": "NONE",
        "features": "NONE",
        "target": "SARSCOV2",
        "match_ref": False,
        "segmented": False,
        "min_coverage": 30,
        "primer_mismatch_rate": 0.1,
        "amplicon_type": "end-to-end",
        "fragment_lookaround_size": None,
        "disable_presets": False,
        "presets": True,
    }
    defaults.update(overrides)
    return Namespace(**defaults)


def test_samplesheet_enforce_absolute_paths_non_dataframe_raises() -> None:
    """Test that non-DataFrame input raises TypeError.

    Verifies that passing a non-DataFrame object (string) to
    samplesheet_enforce_absolute_paths raises TypeError with a descriptive
    message about pandas DataFrame.
    """
    with pytest.raises(TypeError, match="pandas DataFrame"):
        samplesheet_enforce_absolute_paths("not-a-dataframe")  # type: ignore[arg-type]


def test_samplesheet_enforce_absolute_paths_converts_path_columns(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Test conversion of relative paths to absolute paths in DataFrame columns.

    Verifies that samplesheet_enforce_absolute_paths converts relative paths
    to absolute paths in PRIMERS, FEATURES, and REFERENCE columns, handles
    None/NaN values correctly, and preserves non-path columns.

    Parameters
    ----------
    tmp_path : Path
        Temporary directory provided by pytest.
    monkeypatch : pytest.MonkeyPatch
        Pytest fixture for modifying behavior.
    """
    monkeypatch.chdir(tmp_path)
    df = pd.DataFrame(
        {
            "PRIMERS": ["rel/primers.bed", None],
            "FEATURES": ["features.gff", "feature2.gff"],
            "REFERENCE": ["ref.fasta", "ref2.fasta"],
            "VIRUS": ["x", "y"],
        }
    )

    result = samplesheet_enforce_absolute_paths(df)

    assert result.loc[0, "PRIMERS"] == str((tmp_path / "rel/primers.bed").resolve())
    assert pd.isna(result.loc[1, "PRIMERS"])
    assert result.loc[0, "FEATURES"] == str((tmp_path / "features.gff").resolve())
    assert result.loc[0, "REFERENCE"] == str((tmp_path / "ref.fasta").resolve())


def test_file_exists_handles_none_and_real_file(tmp_path: Path) -> None:
    """Test file_exists returns True for NONE and existing files, False otherwise.

    Verifies that the file_exists helper:
    - Returns True for the "NONE" sentinel string
    - Returns True for actual existing file paths
    - Returns False for non-existent paths

    Parameters
    ----------
    tmp_path : Path
        Temporary directory provided by pytest.
    """
    fpath = tmp_path / "a.txt"
    fpath.write_text("x", encoding="utf-8")

    assert file_exists("NONE") is True
    assert file_exists(str(fpath)) is True
    assert file_exists(str(tmp_path / "missing.txt")) is False


@pytest.mark.parametrize(
    "ext,expected_excel,expected_csv,expected_tsv",
    [
        (".xlsx", True, False, False),
        (".xls", True, False, False),
        (".csv", False, True, False),
        (".tsv", False, False, True),
        (".txt", False, False, False),
    ],
)
def test_filetype_helpers(ext: str, expected_excel: bool, expected_csv: bool, expected_tsv: bool) -> None:
    """Test file type detection for Excel, CSV, and TSV extensions.

    Verifies that is_excel_file, is_csv_file, and is_tsv_file correctly
    identify file formats by extension.

    Parameters
    ----------
    ext : str
        File extension to test (e.g., ".xlsx", ".csv", ".tsv").
    expected_excel : bool
        Expected result for is_excel_file check.
    expected_csv : bool
        Expected result for is_csv_file check.
    expected_tsv : bool
        Expected result for is_tsv_file check.
    """
    assert is_excel_file(ext) is expected_excel
    assert is_csv_file(ext) is expected_csv
    assert is_tsv_file(ext) is expected_tsv


def test_open_sample_sheet_csv_and_tsv(tmp_path: Path) -> None:
    """Test reading CSV and TSV sample sheet formats.

    Verifies that open_sample_sheet correctly reads both CSV and TSV
    files, preserving column headers and sample data.

    Parameters
    ----------
    tmp_path : Path
        Temporary directory provided by pytest.
    """
    csv_file = tmp_path / "samples.csv"
    csv_file.write_text("SAMPLE,VIRUS,REFERENCE\ns1,v1,ref.fa\n", encoding="utf-8")

    tsv_file = tmp_path / "samples.tsv"
    tsv_file.write_text("SAMPLE\tVIRUS\tREFERENCE\ns2\tv2\tref2.fa\n", encoding="utf-8")

    csv_df = open_sample_sheet(str(csv_file))
    tsv_df = open_sample_sheet(str(tsv_file))

    assert list(csv_df.columns) == ["SAMPLE", "VIRUS", "REFERENCE"]
    assert list(tsv_df.columns) == ["SAMPLE", "VIRUS", "REFERENCE"]
    assert csv_df.iloc[0]["SAMPLE"] == "s1"
    assert tsv_df.iloc[0]["SAMPLE"] == "s2"


def test_open_sample_sheet_empty_file_exits(tmp_path: Path) -> None:
    """Test that empty sample sheet causes SystemExit.

    Verifies that open_sample_sheet exits gracefully when given an
    empty file (no headers or data).

    Parameters
    ----------
    tmp_path : Path
        Temporary directory provided by pytest.

    Raises
    ------
    SystemExit
        When sample sheet file is empty.
    """
    empty = tmp_path / "empty.csv"
    empty.write_text("", encoding="utf-8")

    with pytest.raises(SystemExit):
        open_sample_sheet(str(empty))


def test_open_sample_sheet_reader_error_returns_empty_df(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Test that parsing errors in sample sheet return empty DataFrame.

    Verifies that when pd.read_csv raises an exception, open_sample_sheet
    gracefully returns an empty DataFrame instead of propagating the error.

    Parameters
    ----------
    tmp_path : Path
        Temporary directory provided by pytest.
    monkeypatch : pytest.MonkeyPatch
        Pytest fixture for mocking behavior.
    """
    csv_file = tmp_path / "samples.csv"
    csv_file.write_text("SAMPLE,VIRUS,REFERENCE\ns1,v1,ref.fa\n", encoding="utf-8")

    def _boom(*args, **kwargs):
        raise ValueError("broken parser")

    monkeypatch.setattr("ViroConstrictor.parser.pd.read_csv", _boom)
    result = open_sample_sheet(str(csv_file))

    assert isinstance(result, pd.DataFrame)
    assert result.empty


def test_required_cols_case_insensitive() -> None:
    """Test case-insensitive validation of required columns.

    Verifies that required_cols accepts column names in lowercase and
    validates that all required columns (SAMPLE, VIRUS, REFERENCE) are
    present.
    """
    assert required_cols(["sample", "virus", "reference"]) is True
    assert required_cols(["sample", "virus"]) is False


def test_check_samplesheet_columns_true_and_false() -> None:
    """Test samplesheet column validation for valid and invalid schemas.

    Verifies that check_samplesheet_columns returns True for dataframes
    with required columns and False for incomplete column sets.
    """
    ok = pd.DataFrame(columns=["SAMPLE", "VIRUS", "REFERENCE"])
    bad = pd.DataFrame(columns=["SAMPLE", "VIRUS"])

    assert check_samplesheet_columns(ok) is True
    assert check_samplesheet_columns(bad) is False


def test_check_samplesheet_empty_rows_drops_all_nan_rows() -> None:
    """Test removal of rows with all NaN values from samplesheet.

    Verifies that check_samplesheet_empty_rows filters out rows where
    all values in required columns are NaN/None and preserves valid rows.
    """
    df = pd.DataFrame(
        [
            {"SAMPLE": "s1", "VIRUS": "v1", "REFERENCE": "r1.fa"},
            {"SAMPLE": None, "VIRUS": None, "REFERENCE": None},
        ]
    )

    cleaned = check_samplesheet_empty_rows(df)

    assert len(cleaned) == 1
    assert cleaned.iloc[0]["SAMPLE"] == "s1"


def test_check_samplesheet_rows_valid_minimal(tmp_path: Path) -> None:
    """Test validation of minimal valid samplesheet row.

    Verifies that a minimal samplesheet with required columns (SAMPLE,
    VIRUS, REFERENCE) and valid file paths passes validation.

    Parameters
    ----------
    tmp_path : Path
        Temporary directory provided by pytest.
    """
    ref = tmp_path / "ref.fasta"
    ref.write_text(">ref\nACTG\n", encoding="utf-8")
    df = pd.DataFrame(
        {
            "SAMPLE": ["sample1"],
            "VIRUS": ["virus1"],
            "PRIMERS": ["NONE"],
            "REFERENCE": [str(ref)],
            "FEATURES": ["NONE"],
            "MATCH-REF": [False],
            "PRESET": ["DEFAULT"],
            "PRESET_SCORE": [0.0],
        }
    )

    result = check_samplesheet_rows(df)

    assert not result.empty
    assert result.iloc[0]["SAMPLE"] == "sample1"


def test_check_samplesheet_rows_disallowed_characters_exits(tmp_path: Path) -> None:
    """Test rejection of sample names with disallowed characters.

    Verifies that sample names containing spaces or other invalid characters
    are rejected with a SystemExit.

    Parameters
    ----------
    tmp_path : Path
        Temporary directory provided by pytest.

    Raises
    ------
    SystemExit
        When sample name contains spaces or other invalid characters.
    """
    ref = tmp_path / "ref.fasta"
    ref.write_text(">ref\nACTG\n", encoding="utf-8")
    df = pd.DataFrame(
        {
            "SAMPLE": ["sample bad"],
            "VIRUS": ["virus1"],
            "PRIMERS": ["NONE"],
            "REFERENCE": [str(ref)],
            "FEATURES": ["NONE"],
            "MATCH-REF": [False],
            "PRESET": ["DEFAULT"],
            "PRESET_SCORE": [0.0],
        }
    )

    with pytest.raises(SystemExit):
        check_samplesheet_rows(df)


def test_check_samplesheet_rows_missing_reference_file_exits(tmp_path: Path) -> None:
    """Test rejection of missing reference file.

    Verifies that a samplesheet row with a non-existent reference file
    is rejected with a SystemExit.

    Parameters
    ----------
    tmp_path : Path
        Temporary directory provided by pytest.

    Raises
    ------
    SystemExit
        When reference file does not exist.
    """
    df = pd.DataFrame(
        {
            "SAMPLE": ["sample1"],
            "VIRUS": ["virus1"],
            "PRIMERS": ["NONE"],
            "REFERENCE": [str(tmp_path / "missing.fasta")],
            "FEATURES": ["NONE"],
            "MATCH-REF": [False],
            "PRESET": ["DEFAULT"],
            "PRESET_SCORE": [0.0],
        }
    )

    with pytest.raises(SystemExit):
        check_samplesheet_rows(df)


def test_check_file_extension_valid_invalid_and_none() -> None:
    """Test file extension validation for allowed, invalid, and NONE values.

    Verifies that check_file_extension returns True for:
    - The "NONE" sentinel string
    - Files with allowed extensions
    And returns False for files with disallowed extensions.
    """
    allowed = [".fasta", ".fa"]
    assert check_file_extension(allowed, "NONE") is True
    assert check_file_extension(allowed, "sample.fa") is True
    assert check_file_extension(allowed, "sample.txt") is False


def test_dir_path_and_check_input_files(tmp_path: Path) -> None:
    """Test directory path validation and input file checking.

    Verifies that dir_path validates directory existence and
    CheckInputFiles verifies that a directory contains sequencing files
    in supported formats (.fastq, .fastq.gz, .fq, .fq.gz).

    Parameters
    ----------
    tmp_path : Path
        Temporary directory provided by pytest.
    """
    valid_dir = tmp_path / "input"
    valid_dir.mkdir()
    (valid_dir / "a.fastq.gz").write_text("@id\nA\n+\n!\n", encoding="utf-8")

    invalid_dir = tmp_path / "invalid"
    invalid_dir.mkdir()
    (invalid_dir / "notes.txt").write_text("x", encoding="utf-8")

    assert dir_path(str(valid_dir)) is True
    assert dir_path(str(tmp_path / "missing")) is False
    assert CheckInputFiles(str(valid_dir)) is True
    assert CheckInputFiles(str(invalid_dir)) is False


def test_args_to_df_populates_existing_df_and_presets(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Test conversion of CLI args to DataFrame with preset matching.

    Verifies that args_to_df converts CLI arguments to DataFrame columns,
    resolves file paths to absolute paths, and populates preset information
    from match_preset_name matching.

    Parameters
    ----------
    tmp_path : Path
        Temporary directory provided by pytest.
    monkeypatch : pytest.MonkeyPatch
        Pytest fixture for mocking behavior.
    """
    args = _build_args(tmp_path, primers="NONE", features="NONE")
    existing_df = pd.DataFrame({"INPUTFILE": ["reads.fastq"]}, index=["sample1"])

    monkeypatch.setattr("ViroConstrictor.parser.match_preset_name", lambda target, use_presets: ("DEFAULT", 0.42))

    result = args_to_df(args, existing_df)

    assert result.loc["sample1", "VIRUS"] == "SARSCOV2"
    assert result.loc["sample1", "REFERENCE"] == str(Path(args.reference).resolve())
    assert result.loc["sample1", "PRIMERS"] == "NONE"
    assert result.loc["sample1", "PRESET"] == "DEFAULT"
    assert result.loc["sample1", "PRESET_SCORE"] == pytest.approx(0.42)


def test_sampledir_to_df_illumina_and_nanopore() -> None:
    """Test conversion of sample directories to DataFrames for Illumina and Nanopore.

    Verifies that sampledir_to_df correctly structures data for different
    sequencing platforms:
    - Illumina paired-end (R1, R2 columns)
    - Illumina single-end (INPUTFILE column)
    - Nanopore (INPUTFILE column)
    """
    illumina = {"s1": {"R1": "r1.fastq.gz", "R2": "r2.fastq.gz"}}
    illumina_df = sampledir_to_df(illumina, "illumina")
    assert list(illumina_df.columns) == ["R1", "R2"]

    illumina_single = {"s2": {"R1": "r1.fastq.gz"}}
    illumina_single_df = sampledir_to_df(illumina_single, "illumina")
    assert list(illumina_single_df.columns) == ["INPUTFILE"]

    nanopore = {"s3": "reads.fastq.gz"}
    nanopore_df = sampledir_to_df(nanopore, "nanopore")
    assert list(nanopore_df.columns) == ["INPUTFILE"]


def test_sampledir_to_df_invalid_platform_raises() -> None:
    """Test rejection of unsupported sequencing platform.

    Verifies that sampledir_to_df raises ValueError when given an
    unsupported platform string.

    Raises
    ------
    ValueError
        When platform is not supported.
    """
    with pytest.raises(ValueError, match="not supported"):
        sampledir_to_df({"s1": "a.fastq"}, "pacbio")


def test_convert_log_text_formats_multisample_dict() -> None:
    """Test formatting of multi-sample dictionary as log text.

    Verifies that convert_log_text creates readable log output containing
    sample names and their associated attributes (VIRUS, REFERENCE, etc).
    """
    samples = {
        "sample1": {"VIRUS": "v1", "REFERENCE": "ref1.fa"},
        "sample2": {"VIRUS": "v2", "REFERENCE": "ref2.fa"},
    }

    text = convert_log_text(samples)

    assert "sample1:" in text
    assert "VIRUS=v1" in text
    assert "sample2:" in text
    assert "REFERENCE=ref2.fa" in text


def test_samplesheet_handle_presets_uses_disable_presets_per_row(monkeypatch: pytest.MonkeyPatch) -> None:
    """Test preset handling respects per-row DISABLE-PRESETS column.

    Verifies that _samplesheet_handle_presets reads per-sample DISABLE-PRESETS
    column values and passes the correct use_presets flag to match_preset_name.

    Parameters
    ----------
    monkeypatch : pytest.MonkeyPatch
        Pytest fixture for mocking behavior.
    """
    parser_obj = CLIparser.__new__(CLIparser)
    parser_obj.flags = Namespace(presets=True)

    calls: list[tuple[str, bool]] = []

    def _fake_match(virus: str, use_presets: bool) -> tuple[str, float]:
        calls.append((virus, use_presets))
        return (f"P_{virus}_{use_presets}", 1.0 if use_presets else 0.0)

    monkeypatch.setattr("ViroConstrictor.parser.match_preset_name", _fake_match)

    df = pd.DataFrame(
        {
            "VIRUS": ["A", "B", "C", "D"],
            "DISABLE-PRESETS": ["TRUE", "FALSE", "maybe", None],
        }
    )

    result = parser_obj._samplesheet_handle_presets(df)

    assert result.loc[0, "DISABLE-PRESETS"] is True
    assert result.loc[1, "DISABLE-PRESETS"] is False
    assert calls == [("A", False), ("B", True), ("C", True), ("D", True)]
    assert result.loc[0, "PRESET"] == "P_A_False"
    assert result.loc[1, "PRESET"] == "P_B_True"


def test_print_missing_asset_warning_no_sheet_exits_when_assets_missing(tmp_path: Path) -> None:
    """Test that missing assets with no samplesheet causes SystemExit.

    Verifies that _print_missing_asset_warning exits when required
    reference/primer/feature files are None and no samplesheet is provided.

    Parameters
    ----------
    tmp_path : Path
        Temporary directory provided by pytest.

    Raises
    ------
    SystemExit
        When required assets are missing and no samplesheet is provided.
    """
    parser_obj = CLIparser.__new__(CLIparser)
    args = _build_args(tmp_path, primers=None, reference=None, features=None, target=None)

    with pytest.raises(SystemExit):
        parser_obj._print_missing_asset_warning(args, sheet_present=False)


def test_print_missing_asset_warning_sheet_present_does_not_exit(tmp_path: Path) -> None:
    """Test that missing assets with samplesheet does not cause exit.

    Verifies that _print_missing_asset_warning does not exit when a
    samplesheet is provided, even with missing assets (allows samplesheet
    to supply defaults).

    Parameters
    ----------
    tmp_path : Path
        Temporary directory provided by pytest.
    """
    parser_obj = CLIparser.__new__(CLIparser)
    args = _build_args(tmp_path)

    parser_obj._print_missing_asset_warning(args, sheet_present=True)


def test_make_samples_dict_exits_when_no_valid_input_files(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Test that missing input files causes SystemExit.

    Verifies that _make_samples_dict exits when the input directory
    contains no valid sequencing files.

    Parameters
    ----------
    tmp_path : Path
        Temporary directory provided by pytest.
    monkeypatch : pytest.MonkeyPatch
        Pytest fixture for mocking behavior.

    Raises
    ------
    SystemExit
        When input directory contains no valid sequence files.
    """
    parser_obj = CLIparser.__new__(CLIparser)
    parser_obj.flags = Namespace(presets=True)
    args = _build_args(tmp_path)

    monkeypatch.setattr("ViroConstrictor.parser.CheckInputFiles", lambda _indir: False)

    with pytest.raises(SystemExit):
        parser_obj._make_samples_dict(None, args, {"s1": "reads.fastq.gz"})


def test_make_samples_dict_without_samplesheet_applies_cli_defaults(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Test CLI defaults applied to samples without samplesheet.

    Verifies that when no samplesheet is provided, CLI arguments supply
    default values for VIRUS, REFERENCE, PRIMERS, FEATURES, and PRESET.

    Parameters
    ----------
    tmp_path : Path
        Temporary directory provided by pytest.
    monkeypatch : pytest.MonkeyPatch
        Pytest fixture for mocking behavior.
    """
    parser_obj = CLIparser.__new__(CLIparser)
    parser_obj.flags = Namespace(presets=True)
    args = _build_args(tmp_path)

    monkeypatch.setattr("ViroConstrictor.parser.CheckInputFiles", lambda _indir: True)
    monkeypatch.setattr("ViroConstrictor.parser.match_preset_name", lambda _v, _u: ("DEFAULT", 0.5))

    result = parser_obj._make_samples_dict(None, args, {"sample1": "reads.fastq.gz"})

    assert "sample1" in result
    assert result["sample1"]["VIRUS"] == args.target
    assert result["sample1"]["REFERENCE"] == str(Path(args.reference).resolve())
    assert result["sample1"]["PRIMERS"] == "NONE"
    assert result["sample1"]["FEATURES"] == "NONE"
    assert result["sample1"]["PRESET"] == "DEFAULT"


def test_make_samples_dict_with_samplesheet_fills_defaults(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Test defaults filled from CLI for samples with samplesheet.

    Verifies that missing fields in samplesheet rows are filled with CLI
    defaults (min_coverage, primer_mismatch_rate, amplicon_type, etc).

    Parameters
    ----------
    tmp_path : Path
        Temporary directory provided by pytest.
    monkeypatch : pytest.MonkeyPatch
        Pytest fixture for mocking behavior.
    """
    parser_obj = CLIparser.__new__(CLIparser)
    parser_obj.flags = Namespace(presets=True)

    input_dir = tmp_path / "input"
    input_dir.mkdir()
    args = _build_args(
        tmp_path,
        input=str(input_dir),
        min_coverage=55,
        primer_mismatch_rate=0.2,
        disable_presets=False,
        fragment_lookaround_size=10,
    )

    monkeypatch.setattr("ViroConstrictor.parser.CheckInputFiles", lambda _indir: True)
    monkeypatch.setattr("ViroConstrictor.parser.match_preset_name", lambda _v, use_presets: ("DEFAULT", 1.0 if use_presets else 0.0))
    monkeypatch.setattr("ViroConstrictor.parser.GenBank.is_genbank", lambda _path: False)

    df = pd.DataFrame(
        {
            "SAMPLE": ["sample1"],
            "VIRUS": ["virus1"],
            "REFERENCE": [args.reference],
            "PRIMERS": ["NONE"],
            "FEATURES": ["NONE"],
        }
    )
    filedict = {"sample1": "reads.fastq.gz", "sample2": "reads2.fastq.gz"}

    result = parser_obj._make_samples_dict(df, args, filedict)

    assert "sample1" in result
    assert result["sample1"]["MIN-COVERAGE"] == 55
    assert result["sample1"]["PRIMER-MISMATCH-RATE"] == pytest.approx(0.2)
    assert result["sample1"]["PRESET"] == "DEFAULT"
    assert result["sample1"]["FRAGMENT-LOOKAROUND-SIZE"] is None


def test_make_samples_dict_with_samplesheet_no_sample_overlap_exits(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Test that non-overlapping samples between samplesheet and input causes exit.

    Verifies that _make_samples_dict exits when samplesheet contains no
    samples that are found in the input directory.

    Parameters
    ----------
    tmp_path : Path
        Temporary directory provided by pytest.
    monkeypatch : pytest.MonkeyPatch
        Pytest fixture for mocking behavior.

    Raises
    ------
    SystemExit
        When samplesheet contains no samples found in input directory.
    """
    parser_obj = CLIparser.__new__(CLIparser)
    parser_obj.flags = Namespace(presets=True)
    args = _build_args(tmp_path)

    monkeypatch.setattr("ViroConstrictor.parser.CheckInputFiles", lambda _indir: True)

    df = pd.DataFrame(
        {
            "SAMPLE": ["sample_in_sheet"],
            "VIRUS": ["virus1"],
            "REFERENCE": [args.reference],
            "PRIMERS": ["NONE"],
            "FEATURES": ["NONE"],
        }
    )

    with pytest.raises(SystemExit):
        parser_obj._make_samples_dict(df, args, {"sample_in_dir": "reads.fastq.gz"})


@pytest.mark.xfail(
    reason="Unknown samplesheet columns should be rejected with a user-facing validation error; current implementation raises KeyError first"
)
def test_check_samplesheet_rows_unknown_column_is_rejected(tmp_path: Path) -> None:
    """Test rejection of unknown columns in samplesheet.

    Verifies that samplesheet with unrecognized column names should cause
    exit (note: currently fails to exit gracefully, raising KeyError instead).

    Parameters
    ----------
    tmp_path : Path
        Temporary directory provided by pytest.

    Raises
    ------
    SystemExit
        When samplesheet contains unrecognized column names.
    """
    ref = tmp_path / "ref.fasta"
    ref.write_text(">ref\nACTG\n", encoding="utf-8")
    df = pd.DataFrame(
        {
            "SAMPLE": ["sample1"],
            "VIRUS": ["virus1"],
            "PRIMERS": ["NONE"],
            "REFERENCE": [str(ref)],
            "FEATURES": ["NONE"],
            "MATCH-REF": [False],
            "PRESET": ["DEFAULT"],
            "PRESET_SCORE": [0.0],
            "UNKNOWN": ["x"],
        }
    )

    with pytest.raises(SystemExit):
        check_samplesheet_rows(df)


def test_check_samplesheet_rows_required_value_null_exits(tmp_path: Path) -> None:
    """Test rejection of null values in required columns.

    Verifies that rows with null/None values in required columns
    (SAMPLE, VIRUS, REFERENCE) are rejected with SystemExit.

    Parameters
    ----------
    tmp_path : Path
        Temporary directory provided by pytest.

    Raises
    ------
    SystemExit
        When required column contains null value.
    """
    ref = tmp_path / "ref.fasta"
    ref.write_text(">ref\nACTG\n", encoding="utf-8")
    df = pd.DataFrame(
        {
            "SAMPLE": [None],
            "VIRUS": ["virus1"],
            "PRIMERS": ["NONE"],
            "REFERENCE": [str(ref)],
            "FEATURES": ["NONE"],
            "MATCH-REF": [False],
            "PRESET": ["DEFAULT"],
            "PRESET_SCORE": [0.0],
        }
    )

    with pytest.raises(SystemExit):
        check_samplesheet_rows(df)


def test_check_samplesheet_rows_required_dtype_mismatch_exits(tmp_path: Path) -> None:
    """Test rejection of incorrect data type in required columns.

    Verifies that required columns with incorrect data types (e.g., int
    instead of string for SAMPLE) are rejected with SystemExit.

    Parameters
    ----------
    tmp_path : Path
        Temporary directory provided by pytest.

    Raises
    ------
    SystemExit
        When required column has wrong data type.
    """
    ref = tmp_path / "ref.fasta"
    ref.write_text(">ref\nACTG\n", encoding="utf-8")
    df = pd.DataFrame(
        {
            "SAMPLE": [1],
            "VIRUS": ["virus1"],
            "PRIMERS": ["NONE"],
            "REFERENCE": [str(ref)],
            "FEATURES": ["NONE"],
            "MATCH-REF": [False],
            "PRESET": ["DEFAULT"],
            "PRESET_SCORE": [0.0],
        }
    )

    with pytest.raises(SystemExit):
        check_samplesheet_rows(df)


@pytest.mark.parametrize(
    "filename,allowed,expected",
    [
        ("reads.fastq.gz", [".fastq", ".fastq.gz"], True),
        ("reads.fq.gz", [".fastq", ".fastq.gz"], False),
        ("reads.fa", [".fasta", ".fa"], True),
        ("reads.fasta.gz", [".fasta", ".fa"], False),
    ],
)
def test_check_file_extension_edge_cases(filename: str, allowed: list[str], expected: bool) -> None:
    """Test file extension validation with compound and edge case extensions.

    Parameters
    ----------
    filename : str
        Filename to validate.
    allowed : list[str]
        List of allowed extensions.
    expected : bool
        Expected validation result.
    """
    assert check_file_extension(allowed, filename) is expected


def test_samplesheet_handle_presets_without_disable_column_uses_cli_flag(monkeypatch: pytest.MonkeyPatch) -> None:
    """Test preset handling uses CLI flag when DISABLE-PRESETS column is absent.

    Verifies that when the samplesheet lacks DISABLE-PRESETS column,
    the parser uses the CLI flag (flags.presets) for all samples.

    Parameters
    ----------
    monkeypatch : pytest.MonkeyPatch
        Pytest fixture for mocking behavior.
    """
    parser_obj = CLIparser.__new__(CLIparser)
    parser_obj.flags = Namespace(presets=False)

    calls: list[tuple[str, bool]] = []

    def _fake_match(virus: str, use_presets: bool) -> tuple[str, float]:
        calls.append((virus, use_presets))
        return ("DEFAULT", 0.0)

    monkeypatch.setattr("ViroConstrictor.parser.match_preset_name", _fake_match)

    df = pd.DataFrame({"VIRUS": ["A", "B"]})
    result = parser_obj._samplesheet_handle_presets(df)

    assert list(result["PRESET"]) == ["DEFAULT", "DEFAULT"]
    assert calls == [("A", False), ("B", False)]


@pytest.mark.parametrize(
    "input_exists,reference,reference_exists,reference_ext,primers,primers_exists,primers_ext,features,features_exists,features_ext,scheduler_valid,expected_fragments",
    [
        (False, "ref.fa", True, True, "p.fa", True, True, "f.gff", True, True, True, ["not a directory"]),
        (True, "NONE", True, True, "p.fa", True, True, "f.gff", True, True, True, ["cannot be given for the reference"]),
        (
            True,
            "ref.txt",
            True,
            False,
            "p.txt",
            True,
            False,
            "f.txt",
            True,
            False,
            True,
            ["reference", "primers", "features"],
        ),
        (
            True,
            "ref.fa",
            False,
            True,
            "p.fa",
            False,
            True,
            "f.gff",
            False,
            True,
            False,
            ["not an existing file", "not a valid scheduler"],
        ),
    ],
)
def test_validate_cli_args_collects_unhappy_errors(
    monkeypatch: pytest.MonkeyPatch,
    input_exists: bool,
    reference: str,
    reference_exists: bool,
    reference_ext: bool,
    primers: str,
    primers_exists: bool,
    primers_ext: bool,
    features: str,
    features_exists: bool,
    features_ext: bool,
    scheduler_valid: bool,
    expected_fragments: list[str],
) -> None:
    """Test CLI argument validation collects multiple errors.

    Parameters
    ----------
    input_exists : bool
        Whether input directory exists.
    reference : str
        Reference file path.
    reference_exists : bool
        Whether reference file exists.
    reference_ext : bool
        Whether reference file has valid extension.
    primers : str
        Primers file path.
    primers_exists : bool
        Whether primers file exists.
    primers_ext : bool
        Whether primers file has valid extension.
    features : str
        Features file path.
    features_exists : bool
        Whether features file exists.
    features_ext : bool
        Whether features file has valid extension.
    scheduler_valid : bool
        Whether scheduler configuration is valid.
    expected_fragments : list[str]
        Expected error message fragments.
    """
    parser_obj = CLIparser.__new__(CLIparser)
    parser_obj.flags = Namespace(
        input="/tmp/in",
        samplesheet="sheet.csv",
        reference=reference,
        primers=primers,
        features=features,
        scheduler="invalid",
    )

    monkeypatch.setattr("ViroConstrictor.parser.dir_path", lambda _p: input_exists)

    def _file_exists(path: str) -> bool:
        mapping = {
            "sheet.csv": True,
            reference: reference_exists,
            primers: primers_exists,
            features: features_exists,
        }
        return mapping.get(path, False)

    monkeypatch.setattr("ViroConstrictor.parser.file_exists", _file_exists)

    def _ext_check(allowed_extensions: list[str], fname: str) -> bool:
        if fname == reference:
            return reference_ext
        if fname == primers:
            return primers_ext
        if fname == features:
            return features_ext
        return True

    monkeypatch.setattr("ViroConstrictor.parser.check_file_extension", _ext_check)
    monkeypatch.setattr("ViroConstrictor.parser.Scheduler.is_valid", lambda _sched: scheduler_valid)
    monkeypatch.setattr("ViroConstrictor.parser.Scheduler.supported_schedulers", lambda: ["none", "slurm"])

    errors = parser_obj._validate_cli_args()

    assert errors is not None
    combined_errors = "\n".join(errors)
    for fragment in expected_fragments:
        assert fragment in combined_errors


def test_check_sample_sheet_missing_required_columns_exits(monkeypatch: pytest.MonkeyPatch) -> None:
    """Test that missing required samplesheet columns causes exit.

    Raises
    ------
    SystemExit
        When samplesheet lacks required columns.
    """
    parser_obj = CLIparser.__new__(CLIparser)
    parser_obj.flags = Namespace(presets=True)

    monkeypatch.setattr("ViroConstrictor.parser.open_sample_sheet", lambda _f: pd.DataFrame({"SAMPLE": ["s1"]}))
    monkeypatch.setattr("ViroConstrictor.parser.check_samplesheet_columns", lambda _df: False)

    with pytest.raises(SystemExit):
        parser_obj._check_sample_sheet("dummy.csv")


def test_check_sample_sheet_empty_returns_empty_df(monkeypatch: pytest.MonkeyPatch) -> None:
    """Test that empty samplesheet returns empty DataFrame.

    Verifies that _check_sample_sheet returns an empty DataFrame when
    the input file contains no rows.

    Parameters
    ----------
    monkeypatch : pytest.MonkeyPatch
        Pytest fixture for mocking behavior.
    """
    parser_obj = CLIparser.__new__(CLIparser)
    parser_obj.flags = Namespace(presets=True)
    monkeypatch.setattr("ViroConstrictor.parser.open_sample_sheet", lambda _f: pd.DataFrame())

    result = parser_obj._check_sample_sheet("dummy.csv")

    assert result.empty


def test_check_sample_properties_missing_reference_exits() -> None:
    """Test that missing sample reference file causes exit.

    Raises
    ------
    SystemExit
        When sample reference file does not exist.
    """
    parser_obj = CLIparser.__new__(CLIparser)
    sampleinfo = {"s1": {"REFERENCE": str(Path(tempfile.gettempdir()) / f"missing_{uuid.uuid4().hex}.fasta")}}

    with pytest.raises(SystemExit):
        parser_obj._check_sample_properties(sampleinfo)


def test_check_sample_properties_checks_unique_references(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Test that unique reference files are validated once.

    Verifies that when multiple samples reference the same file,
    CheckReferenceFile is called only once for that unique reference.

    Parameters
    ----------
    tmp_path : Path
        Temporary directory provided by pytest.
    monkeypatch : pytest.MonkeyPatch
        Pytest fixture for mocking behavior.
    """
    parser_obj = CLIparser.__new__(CLIparser)
    ref = tmp_path / "ref.fasta"
    ref.write_text(">r\nACTG\n", encoding="utf-8")
    sampleinfo = {
        "s1": {"REFERENCE": str(ref)},
        "s2": {"REFERENCE": str(ref)},
    }
    calls: list[str] = []

    monkeypatch.setattr("ViroConstrictor.parser.CheckReferenceFile", lambda p: calls.append(p))

    parser_obj._check_sample_properties(sampleinfo)

    assert calls == [str(ref)]


def test_make_samples_dict_virus_empty_exits(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Test that empty virus name causes exit.

    Raises
    ------
    SystemExit
        When virus field is empty string.
    """
    parser_obj = CLIparser.__new__(CLIparser)
    parser_obj.flags = Namespace(presets=True)
    args = _build_args(tmp_path)

    monkeypatch.setattr("ViroConstrictor.parser.CheckInputFiles", lambda _indir: True)

    df = pd.DataFrame(
        {
            "SAMPLE": ["sample1"],
            "VIRUS": [""],
            "REFERENCE": [args.reference],
            "PRIMERS": ["NONE"],
            "FEATURES": ["NONE"],
        }
    )

    with pytest.raises(SystemExit):
        parser_obj._make_samples_dict(df, args, {"sample1": "reads.fastq.gz"})


def test_make_samples_dict_missing_default_reference_exits(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Test that missing default reference file causes exit.

    Verifies that _make_samples_dict exits when reference is None (not
    provided via CLI) and not supplied in samplesheet row.

    Parameters
    ----------
    tmp_path : Path
        Temporary directory provided by pytest.
    monkeypatch : pytest.MonkeyPatch
        Pytest fixture for mocking behavior.

    Raises
    ------
    SystemExit
        When reference is None and not provided in samplesheet.
    """
    parser_obj = CLIparser.__new__(CLIparser)
    parser_obj.flags = Namespace(presets=True)
    args = _build_args(tmp_path, reference=None)

    monkeypatch.setattr("ViroConstrictor.parser.CheckInputFiles", lambda _indir: True)

    df = pd.DataFrame(
        {
            "SAMPLE": ["sample1"],
            "VIRUS": ["virus1"],
            "REFERENCE": [None],
            "PRIMERS": ["NONE"],
            "FEATURES": ["NONE"],
        }
    )

    with pytest.raises(SystemExit):
        parser_obj._make_samples_dict(df, args, {"sample1": "reads.fastq.gz"})


def test_make_samples_dict_missing_default_primers_exits(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Test that missing default primers file causes exit.

    Verifies that _make_samples_dict exits when primers is None (not
    provided via CLI) and not supplied in samplesheet row.

    Parameters
    ----------
    tmp_path : Path
        Temporary directory provided by pytest.
    monkeypatch : pytest.MonkeyPatch
        Pytest fixture for mocking behavior.

    Raises
    ------
    SystemExit
        When primers is None and not provided in samplesheet.
    """
    parser_obj = CLIparser.__new__(CLIparser)
    parser_obj.flags = Namespace(presets=True)
    args = _build_args(tmp_path, primers=None)

    monkeypatch.setattr("ViroConstrictor.parser.CheckInputFiles", lambda _indir: True)
    monkeypatch.setattr("ViroConstrictor.parser.GenBank.is_genbank", lambda _path: False)

    df = pd.DataFrame(
        {
            "SAMPLE": ["sample1"],
            "VIRUS": ["virus1"],
            "REFERENCE": [args.reference],
            "PRIMERS": [None],
            "FEATURES": ["NONE"],
        }
    )

    with pytest.raises(SystemExit):
        parser_obj._make_samples_dict(df, args, {"sample1": "reads.fastq.gz"})


def test_make_samples_dict_missing_default_features_exits(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Test that missing default features file causes exit.

    Verifies that _make_samples_dict exits when features is None (not
    provided via CLI) and not supplied in samplesheet row.

    Parameters
    ----------
    tmp_path : Path
        Temporary directory provided by pytest.
    monkeypatch : pytest.MonkeyPatch
        Pytest fixture for mocking behavior.

    Raises
    ------
    SystemExit
        When features is None and not provided in samplesheet.
    """
    parser_obj = CLIparser.__new__(CLIparser)
    parser_obj.flags = Namespace(presets=True)
    args = _build_args(tmp_path, features=None)

    monkeypatch.setattr("ViroConstrictor.parser.CheckInputFiles", lambda _indir: True)
    monkeypatch.setattr("ViroConstrictor.parser.GenBank.is_genbank", lambda _path: False)

    df = pd.DataFrame(
        {
            "SAMPLE": ["sample1"],
            "VIRUS": ["virus1"],
            "REFERENCE": [args.reference],
            "PRIMERS": ["NONE"],
            "FEATURES": [None],
        }
    )

    with pytest.raises(SystemExit):
        parser_obj._make_samples_dict(df, args, {"sample1": "reads.fastq.gz"})


def test_make_samples_dict_more_samples_in_sheet_than_input_exits(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Test that duplicate sample names in samplesheet causes exit.

    Verifies that _make_samples_dict exits when samplesheet contains
    duplicate sample names (multiple rows with same SAMPLE value).

    Parameters
    ----------
    tmp_path : Path
        Temporary directory provided by pytest.
    monkeypatch : pytest.MonkeyPatch
        Pytest fixture for mocking behavior.

    Raises
    ------
    SystemExit
        When duplicate sample names exist in samplesheet.
    """
    parser_obj = CLIparser.__new__(CLIparser)
    parser_obj.flags = Namespace(presets=True)
    args = _build_args(tmp_path)

    monkeypatch.setattr("ViroConstrictor.parser.CheckInputFiles", lambda _indir: True)

    df = pd.DataFrame(
        {
            "SAMPLE": ["sample1", "sample1"],
            "VIRUS": ["virus1", "virus2"],
            "REFERENCE": [args.reference, args.reference],
            "PRIMERS": ["NONE", "NONE"],
            "FEATURES": ["NONE", "NONE"],
        }
    )

    with pytest.raises(SystemExit):
        parser_obj._make_samples_dict(df, args, {"sample1": "reads.fastq.gz"})


def test_make_samples_dict_handles_genbank_and_applies_user_intent_defaults(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Test GenBank handling and application of user-specified defaults.

    Verifies that _make_samples_dict correctly:
    - Processes GenBank files (.gbk) referenced in REFERENCE column
    - Applies CLI defaults for missing values (primers, features, presets)
    - Respects user-specified DISABLE-PRESETS and FRAGMENT-LOOKAROUND-SIZE

    Parameters
    ----------
    tmp_path : Path
        Temporary directory provided by pytest.
    monkeypatch : pytest.MonkeyPatch
        Pytest fixture for mocking behavior.
    """
    parser_obj = CLIparser.__new__(CLIparser)
    parser_obj.flags = Namespace(presets=False)

    input_dir = tmp_path / "input"
    input_dir.mkdir()
    args = _build_args(
        tmp_path,
        input=str(input_dir),
        reference=str(tmp_path / "default_ref.fasta"),
        primers="NONE",
        features="NONE",
        match_ref=False,
        segmented=False,
        amplicon_type="end-to-end",
        fragment_lookaround_size=None,
        disable_presets=True,
    )

    df = pd.DataFrame(
        {
            "SAMPLE": ["sample1"],
            "VIRUS": ["virus1"],
            "REFERENCE": [str(tmp_path / "sample_ref.gbk")],
            "PRIMERS": [""],
            "FEATURES": [None],
            "PRIMER-MISMATCH-RATE": [float("nan")],
            "PRESET": [""],
            "PRESET_SCORE": [None],
            "FRAGMENT-LOOKAROUND-SIZE": [10],
            "DISABLE-PRESETS": ["MAYBE"],
        }
    )

    monkeypatch.setattr("ViroConstrictor.parser.CheckInputFiles", lambda _indir: True)
    monkeypatch.setattr(
        "ViroConstrictor.parser.GenBank.is_genbank",
        lambda p: str(p).endswith(".gbk"),
    )
    monkeypatch.setattr(
        "ViroConstrictor.parser.GenBank.split_genbank",
        lambda _path, emit_target: (str(tmp_path / "split_ref.fasta"), str(tmp_path / "split_features.gff"), "virus1"),
    )
    preset_calls: list[bool] = []

    def _fake_match_preset_name(_virus: str, use_presets: bool) -> tuple[str, float]:
        preset_calls.append(use_presets)
        return ("DEFAULT", 0.75 if use_presets else 0.0)

    monkeypatch.setattr("ViroConstrictor.parser.match_preset_name", _fake_match_preset_name)

    result = parser_obj._make_samples_dict(df, args, {"sample1": "reads.fastq.gz"})

    assert result["sample1"]["REFERENCE"] == str(tmp_path / "split_ref.fasta")
    assert result["sample1"]["FEATURES"] == str(tmp_path / "split_features.gff")
    assert result["sample1"]["PRIMERS"] == "NONE"
    assert result["sample1"]["PRIMER-MISMATCH-RATE"] == pytest.approx(0.1)
    assert result["sample1"]["MATCH-REF"] is False
    assert result["sample1"]["SEGMENTED"] is False
    assert result["sample1"]["PRESET"] == "DEFAULT"
    assert result["sample1"]["PRESET_SCORE"] == pytest.approx(0.0)
    assert result["sample1"]["FRAGMENT-LOOKAROUND-SIZE"] is None
    # Invalid DISABLE-PRESETS values should fall back to CLI behavior.
    assert preset_calls
    assert preset_calls[0] is False


def test_make_samples_dict_fragmented_with_invalid_lookaround_uses_default(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Test that invalid fragment lookaround value defaults to CLI value.

    Verifies that invalid FRAGMENT-LOOKAROUND-SIZE values (non-numeric, negative)
    are replaced with the CLI-specified default.

    Parameters
    ----------
    tmp_path : Path
        Temporary directory provided by pytest.
    monkeypatch : pytest.MonkeyPatch
        Pytest fixture for mocking behavior.
    """
    parser_obj = CLIparser.__new__(CLIparser)
    parser_obj.flags = Namespace(presets=True)

    input_dir = tmp_path / "input"
    input_dir.mkdir()
    args = _build_args(
        tmp_path,
        input=str(input_dir),
        amplicon_type="fragmented",
        fragment_lookaround_size=25,
        disable_presets=False,
    )

    preset_calls: list[bool] = []

    df = pd.DataFrame(
        {
            "SAMPLE": ["sample1"],
            "VIRUS": ["virus1"],
            "REFERENCE": [args.reference],
            "PRIMERS": ["NONE"],
            "FEATURES": ["NONE"],
            "FRAGMENT-LOOKAROUND-SIZE": ["not-an-int"],
            "DISABLE-PRESETS": ["FALSE"],
            "PRESET": [""],
            "PRESET_SCORE": [None],
        }
    )

    monkeypatch.setattr("ViroConstrictor.parser.CheckInputFiles", lambda _indir: True)
    monkeypatch.setattr("ViroConstrictor.parser.GenBank.is_genbank", lambda _path: False)

    def _fake_match_preset_name(_virus: str, use_presets: bool) -> tuple[str, float]:
        preset_calls.append(use_presets)
        return ("DEFAULT", 1.0 if use_presets else 0.0)

    monkeypatch.setattr("ViroConstrictor.parser.match_preset_name", _fake_match_preset_name)

    result = parser_obj._make_samples_dict(df, args, {"sample1": "reads.fastq.gz"})

    assert result["sample1"]["FRAGMENT-LOOKAROUND-SIZE"] == 25
    # Explicit FALSE means presets are enabled for matching.
    assert preset_calls
    assert preset_calls[0] is True


def test_get_args_with_no_arguments_exits() -> None:
    """Test that parser exits when no arguments provided.

    Raises
    ------
    SystemExit
        When called with empty argument list.
    """
    parser_obj = CLIparser.__new__(CLIparser)

    with pytest.raises(SystemExit):
        parser_obj._get_args([])


def test_get_args_parses_required_values(tmp_path: Path) -> None:
    """Test parsing of required CLI arguments.

    Verifies that _get_args correctly parses input, output, platform,
    amplicon-type, and other required CLI arguments into a namespace.

    Parameters
    ----------
    tmp_path : Path
        Temporary directory provided by pytest.
    """
    parser_obj = CLIparser.__new__(CLIparser)
    out_dir = tmp_path / "out"

    flags = parser_obj._get_args(
        [
            "--input",
            str(tmp_path),
            "--output",
            str(out_dir),
            "--platform",
            "nanopore",
            "--amplicon-type",
            "end-to-end",
        ]
    )

    assert flags.input == str(tmp_path)
    assert flags.output == str(out_dir)
    assert flags.platform == "nanopore"
    assert flags.amplicon_type == "end-to-end"
    assert flags.scheduler == "auto"
    assert flags.min_coverage == 30
    assert flags.disable_presets is False


def test_parse_genbank_updates_flags(monkeypatch: pytest.MonkeyPatch) -> None:
    """Test that GenBank file parsing updates parser flags.

    Verifies that parse_genbank extracts references, features, and target
    from a GenBank file and updates parser flags with the results.

    Parameters
    ----------
    monkeypatch : pytest.MonkeyPatch
        Pytest fixture for mocking behavior.
    """
    parser_obj = CLIparser.__new__(CLIparser)
    parser_obj.flags = Namespace(reference="old-ref", features="old-features", target="old-target")

    monkeypatch.setattr(
        "ViroConstrictor.parser.GenBank.split_genbank",
        lambda _path, emit_target: ("new-ref.fasta", "new-features.gff", "new-target"),
    )

    parser_obj.parse_genbank("input.gbk")

    assert parser_obj.flags.reference == "new-ref.fasta"
    assert parser_obj.flags.features == "new-features.gff"
    assert parser_obj.flags.target == "new-target"


def test_validate_cli_args_samplesheet_missing_and_invalid_extension(monkeypatch: pytest.MonkeyPatch) -> None:
    """Test detection of missing samplesheet with invalid file extension.

    Verifies that _validate_cli_args collects errors for:
    - Samplesheet file not existing
    - Samplesheet file having invalid extension

    Parameters
    ----------
    monkeypatch : pytest.MonkeyPatch
        Pytest fixture for mocking behavior.
    """
    parser_obj = CLIparser.__new__(CLIparser)
    parser_obj.flags = Namespace(
        input="/tmp/input",
        samplesheet="sheet.bad",
        reference=None,
        primers=None,
        features=None,
        scheduler="auto",
    )

    monkeypatch.setattr("ViroConstrictor.parser.dir_path", lambda _p: True)
    monkeypatch.setattr("ViroConstrictor.parser.file_exists", lambda _p: False)
    monkeypatch.setattr("ViroConstrictor.parser.check_file_extension", lambda **_kwargs: False)
    monkeypatch.setattr("ViroConstrictor.parser.Scheduler.is_valid", lambda _sched: True)

    errors = parser_obj._validate_cli_args()

    assert errors is not None
    assert any("not an existing file" in msg for msg in errors)
    assert any("valid file extension" in msg for msg in errors)


def test_get_paths_for_workflow_creates_output_directory(tmp_path: Path) -> None:
    """Test workflow path setup creates output directory.

    Verifies that _get_paths_for_workflow:
    - Creates output directory if it doesn't exist
    - Returns absolute paths for input, workdir, exec_start
    - Returns correct snakefile paths for main and match-ref workflows

    Parameters
    ----------
    tmp_path : Path
        Temporary directory provided by pytest.
    """
    parser_obj = CLIparser.__new__(CLIparser)
    output_dir = tmp_path / "new-output"
    flags = Namespace(input=str(tmp_path), output=str(output_dir))

    input_path, workdir, exec_start, snakefile, match_ref_snakefile = parser_obj._get_paths_for_workflow(flags)

    assert Path(input_path).is_absolute()
    assert Path(workdir) == output_dir.resolve()
    assert output_dir.exists()
    assert Path(exec_start).is_absolute()
    assert snakefile.endswith("workflow/main/workflow.smk")
    assert match_ref_snakefile.endswith("workflow/match_ref/workflow.smk")


def test_check_sample_sheet_happy_path_returns_processed_dataframe(monkeypatch: pytest.MonkeyPatch) -> None:
    """Test successful samplesheet processing through all validation steps.

    Verifies that _check_sample_sheet:
    - Opens and validates the samplesheet
    - Checks columns, enforces absolute paths, validates rows
    - Adds preset information when presets are enabled

    Parameters
    ----------
    monkeypatch : pytest.MonkeyPatch
        Pytest fixture for mocking behavior.
    """
    parser_obj = CLIparser.__new__(CLIparser)
    parser_obj.flags = Namespace(presets=True)

    source_df = pd.DataFrame({"sample": ["s1"], "virus": ["v1"], "reference": ["r.fa"]})
    validated_df = pd.DataFrame({"SAMPLE": ["s1"], "VIRUS": ["v1"], "REFERENCE": ["r.fa"]})

    def _open_sheet(_file: str) -> pd.DataFrame:
        return source_df.copy()

    def _check_columns(df: pd.DataFrame) -> bool:
        assert list(df.columns) == ["SAMPLE", "VIRUS", "REFERENCE"]
        return True

    def _drop_empty(df: pd.DataFrame) -> pd.DataFrame:
        return df

    def _enforce_paths(df: pd.DataFrame) -> pd.DataFrame:
        return df

    def _check_rows(df: pd.DataFrame) -> pd.DataFrame:
        return validated_df

    def _handle_presets(df: pd.DataFrame) -> pd.DataFrame:
        return df.assign(PRESET="DEFAULT", PRESET_SCORE=0.0)

    monkeypatch.setattr("ViroConstrictor.parser.open_sample_sheet", _open_sheet)
    monkeypatch.setattr("ViroConstrictor.parser.check_samplesheet_columns", _check_columns)
    monkeypatch.setattr("ViroConstrictor.parser.check_samplesheet_empty_rows", _drop_empty)
    monkeypatch.setattr("ViroConstrictor.parser.samplesheet_enforce_absolute_paths", _enforce_paths)
    monkeypatch.setattr("ViroConstrictor.parser.check_samplesheet_rows", _check_rows)
    monkeypatch.setattr(parser_obj, "_samplesheet_handle_presets", _handle_presets)

    result = parser_obj._check_sample_sheet("sheet.csv")

    assert result.loc[0, "SAMPLE"] == "s1"
    assert result.loc[0, "PRESET"] == "DEFAULT"
    assert result.loc[0, "PRESET_SCORE"] == pytest.approx(0.0)
    assert {"SAMPLE", "VIRUS", "REFERENCE", "PRESET", "PRESET_SCORE"}.issubset(set(result.columns))


def test_cli_parser_init_with_samplesheet_success(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """Test successful CLIparser initialization with valid samplesheet.

    Verifies that CLIparser.__init__:
    - Parses command-line arguments
    - Validates arguments
    - Loads configuration
    - Initializes scheduler
    - Processes samplesheet
    - Creates samples dictionary

    Parameters
    ----------
    monkeypatch : pytest.MonkeyPatch
        Pytest fixture for mocking behavior.
    tmp_path : Path
        Temporary directory provided by pytest.
    """
    input_dir = tmp_path / "input"
    input_dir.mkdir()
    output_dir = tmp_path / "output"

    flags = Namespace(
        input=str(input_dir),
        output=str(output_dir),
        verbose=False,
        samplesheet="sheet.csv",
        scheduler="auto",
        dryrun=False,
        disable_presets=False,
        platform="nanopore",
        reference="ref.fasta",
        features="NONE",
        target="virus",
    )

    cfg = ConfigParser()
    cfg["COMPUTING"] = {"compmode": "local", "scheduler": "none"}

    monkeypatch.setattr(CLIparser, "_get_args", lambda self, _args: flags)
    monkeypatch.setattr("ViroConstrictor.parser.setup_logger", lambda _out: "logfile")
    monkeypatch.setattr(CLIparser, "_validate_cli_args", lambda self: [])
    monkeypatch.setattr("ViroConstrictor.parser.ReadConfig", lambda _path: cfg)
    monkeypatch.setattr("ViroConstrictor.parser.Scheduler.determine_scheduler", lambda *_args: scheduler_module.Scheduler.LOCAL)
    monkeypatch.setattr(CLIparser, "_print_missing_asset_warning", lambda self, _args, _present: None)
    monkeypatch.setattr(CLIparser, "_check_sample_sheet", lambda self, _sheet: pd.DataFrame({"SAMPLE": ["s1"]}))
    monkeypatch.setattr("ViroConstrictor.parser.GetSamples", lambda _inp, _plat: {"s1": "reads.fastq.gz"})
    monkeypatch.setattr(CLIparser, "_make_samples_dict", lambda self, _df, _flags, _files: {"s1": {"REFERENCE": str(tmp_path / "r.fa")}})
    monkeypatch.setattr(CLIparser, "_get_paths_for_workflow", lambda self, _flags: ("/in", "/out", "/cwd", "main.smk", "mr.smk"))
    monkeypatch.setattr(CLIparser, "_check_sample_properties", lambda self, _samples: None)

    parser_obj = CLIparser(["--ignored"], "~/.config.ini")

    assert parser_obj.flags.presets is True
    assert parser_obj.scheduler == scheduler_module.Scheduler.LOCAL
    assert parser_obj.samples_dict == {"s1": {"REFERENCE": str(tmp_path / "r.fa")}}
    assert "SAMPLE" in parser_obj.samples_df.columns
    assert "REFERENCE" in parser_obj.samples_df.columns
    assert parser_obj.samples_df.loc[0, "SAMPLE"] == "s1"


def test_cli_parser_init_without_samplesheet_triggers_genbank_path(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """Test CLIparser initialization triggers GenBank parsing when reference is GenBank file.

    Verifies that when no samplesheet is provided but reference is a GenBank file,
    CLIparser detects and triggers parse_genbank() to extract reference and feature information.

    Parameters
    ----------
    monkeypatch : pytest.MonkeyPatch
        Pytest fixture for mocking behavior.
    tmp_path : Path
        Temporary directory provided by pytest.
    """
    input_dir = tmp_path / "input"
    input_dir.mkdir()
    output_dir = tmp_path / "output"

    flags = Namespace(
        input=str(input_dir),
        output=str(output_dir),
        verbose=True,
        samplesheet=None,
        scheduler="auto",
        dryrun=False,
        disable_presets=True,
        platform="nanopore",
        reference=str(tmp_path / "ref.gbk"),
        features="NONE",
        target="virus",
    )

    cfg = ConfigParser()
    cfg["COMPUTING"] = {"compmode": "local", "scheduler": "none"}

    parse_called = {"value": False}

    monkeypatch.setattr(CLIparser, "_get_args", lambda self, _args: flags)
    monkeypatch.setattr("ViroConstrictor.parser.setup_logger", lambda _out: "logfile")
    monkeypatch.setattr(CLIparser, "_validate_cli_args", lambda self: [])
    monkeypatch.setattr("ViroConstrictor.parser.ReadConfig", lambda _path: cfg)
    monkeypatch.setattr("ViroConstrictor.parser.Scheduler.determine_scheduler", lambda *_args: scheduler_module.Scheduler.LOCAL)
    monkeypatch.setattr(CLIparser, "_print_missing_asset_warning", lambda self, _args, _present: None)
    monkeypatch.setattr("ViroConstrictor.parser.GenBank.is_genbank", lambda _path: True)
    monkeypatch.setattr(CLIparser, "parse_genbank", lambda self, _ref: parse_called.__setitem__("value", True))
    monkeypatch.setattr("ViroConstrictor.parser.GetSamples", lambda _inp, _plat: {"s1": "reads.fastq.gz"})
    monkeypatch.setattr(CLIparser, "_make_samples_dict", lambda self, _df, _flags, _files: {"s1": {"REFERENCE": str(tmp_path / "r.fa")}})
    monkeypatch.setattr(CLIparser, "_get_paths_for_workflow", lambda self, _flags: ("/in", "/out", "/cwd", "main.smk", "mr.smk"))
    monkeypatch.setattr(CLIparser, "_check_sample_properties", lambda self, _samples: None)

    parser_obj = CLIparser(["--ignored"], "~/.config.ini")

    assert parse_called["value"] is True
    assert parser_obj.flags.presets is False


def test_cli_parser_init_exits_when_cli_errors_present(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """Test CLIparser initialization exits when validation errors are found.

    Raises
    ------
    SystemExit
        When CLI argument validation fails.
    """
    input_dir = tmp_path / "input"
    input_dir.mkdir()

    flags = Namespace(
        input=str(input_dir),
        output=str(tmp_path / "output"),
        verbose=False,
        samplesheet=None,
        scheduler="auto",
        dryrun=False,
        disable_presets=False,
    )

    monkeypatch.setattr(CLIparser, "_get_args", lambda self, _args: flags)
    monkeypatch.setattr("ViroConstrictor.parser.setup_logger", lambda _out: "logfile")
    monkeypatch.setattr(CLIparser, "_validate_cli_args", lambda self: ["bad arg"])

    with pytest.raises(SystemExit):
        CLIparser(["--ignored"], "~/.config.ini")
