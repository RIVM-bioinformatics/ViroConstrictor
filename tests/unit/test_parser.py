"""Unit tests for parser utilities and CLIparser helper methods.

Coverage focus:
- file/path and samplesheet helpers
- dataframe transformation helpers
- key CLIparser internal methods with mocked dependencies
"""

from argparse import Namespace
from pathlib import Path

import pandas as pd
import pytest

from ViroConstrictor.parser import (
    CLIparser,
    CheckInputFiles,
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
    """Create a parser-like Namespace with defaults for internal method tests."""
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
    with pytest.raises(TypeError, match="pandas DataFrame"):
        samplesheet_enforce_absolute_paths("not-a-dataframe")  # type: ignore[arg-type]


def test_samplesheet_enforce_absolute_paths_converts_path_columns(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
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
    assert is_excel_file(ext) is expected_excel
    assert is_csv_file(ext) is expected_csv
    assert is_tsv_file(ext) is expected_tsv


def test_open_sample_sheet_csv_and_tsv(tmp_path: Path) -> None:
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
    empty = tmp_path / "empty.csv"
    empty.write_text("", encoding="utf-8")

    with pytest.raises(SystemExit):
        open_sample_sheet(str(empty))


def test_open_sample_sheet_reader_error_returns_empty_df(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    csv_file = tmp_path / "samples.csv"
    csv_file.write_text("SAMPLE,VIRUS,REFERENCE\ns1,v1,ref.fa\n", encoding="utf-8")

    def _boom(*args, **kwargs):
        raise ValueError("broken parser")

    monkeypatch.setattr("ViroConstrictor.parser.pd.read_csv", _boom)
    result = open_sample_sheet(str(csv_file))

    assert isinstance(result, pd.DataFrame)
    assert result.empty


def test_required_cols_case_insensitive() -> None:
    assert required_cols(["sample", "virus", "reference"]) is True
    assert required_cols(["sample", "virus"]) is False


def test_check_samplesheet_columns_true_and_false() -> None:
    ok = pd.DataFrame(columns=["SAMPLE", "VIRUS", "REFERENCE"])
    bad = pd.DataFrame(columns=["SAMPLE", "VIRUS"])

    assert check_samplesheet_columns(ok) is True
    assert check_samplesheet_columns(bad) is False


def test_check_samplesheet_empty_rows_drops_all_nan_rows() -> None:
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
    allowed = [".fasta", ".fa"]
    assert check_file_extension(allowed, "NONE") is True
    assert check_file_extension(allowed, "sample.fa") is True
    assert check_file_extension(allowed, "sample.txt") is False


def test_dir_path_and_check_input_files(tmp_path: Path) -> None:
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
    args = _build_args(tmp_path, primers="NONE", features="NONE")
    existing_df = pd.DataFrame({"INPUTFILE": ["reads.fastq"]}, index=["sample1"])

    monkeypatch.setattr("ViroConstrictor.parser.match_preset_name", lambda target, use_presets: ("DEFAULT", 0.42))

    result = args_to_df(args, existing_df)

    assert result.loc["sample1", "VIRUS"] == "SARSCOV2"
    assert result.loc["sample1", "REFERENCE"] == str(Path(args.reference).resolve())
    assert result.loc["sample1", "PRIMERS"] == "NONE"
    assert result.loc["sample1", "PRESET"] == "DEFAULT"
    assert result.loc["sample1", "PRESET_SCORE"] == 0.42


def test_sampledir_to_df_illumina_and_nanopore() -> None:
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
    with pytest.raises(ValueError, match="not supported"):
        sampledir_to_df({"s1": "a.fastq"}, "pacbio")


def test_convert_log_text_formats_multisample_dict() -> None:
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
    parser_obj = CLIparser.__new__(CLIparser)
    args = _build_args(tmp_path, primers=None, reference=None, features=None, target=None)

    with pytest.raises(SystemExit):
        parser_obj._print_missing_asset_warning(args, sheet_present=False)


def test_print_missing_asset_warning_sheet_present_does_not_exit(tmp_path: Path) -> None:
    parser_obj = CLIparser.__new__(CLIparser)
    args = _build_args(tmp_path)

    parser_obj._print_missing_asset_warning(args, sheet_present=True)


def test_make_samples_dict_exits_when_no_valid_input_files(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    parser_obj = CLIparser.__new__(CLIparser)
    parser_obj.flags = Namespace(presets=True)
    args = _build_args(tmp_path)

    monkeypatch.setattr("ViroConstrictor.parser.CheckInputFiles", lambda _indir: False)

    with pytest.raises(SystemExit):
        parser_obj._make_samples_dict(None, args, {"s1": "reads.fastq.gz"})


def test_make_samples_dict_without_samplesheet_uses_args_to_df(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    parser_obj = CLIparser.__new__(CLIparser)
    parser_obj.flags = Namespace(presets=True)
    args = _build_args(tmp_path)

    monkeypatch.setattr("ViroConstrictor.parser.CheckInputFiles", lambda _indir: True)

    def _fake_args_to_df(_args: Namespace, existing_df: pd.DataFrame) -> pd.DataFrame:
        existing_df["VIRUS"] = "X"
        existing_df["REFERENCE"] = _args.reference
        existing_df["PRIMERS"] = _args.primers
        existing_df["FEATURES"] = _args.features
        return existing_df

    monkeypatch.setattr("ViroConstrictor.parser.args_to_df", _fake_args_to_df)

    result = parser_obj._make_samples_dict(None, args, {"sample1": "reads.fastq.gz"})

    assert "sample1" in result
    assert result["sample1"]["VIRUS"] == "X"


def test_make_samples_dict_with_samplesheet_fills_defaults(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
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
    assert result["sample1"]["PRIMER-MISMATCH-RATE"] == 0.2
    assert result["sample1"]["PRESET"] == "DEFAULT"
    assert result["sample1"]["FRAGMENT-LOOKAROUND-SIZE"] is None


def test_make_samples_dict_with_samplesheet_no_sample_overlap_exits(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
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


def test_check_samplesheet_rows_unknown_column_raises_keyerror(tmp_path: Path) -> None:
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

    with pytest.raises(KeyError, match="UNKNOWN"):
        check_samplesheet_rows(df)


def test_check_samplesheet_rows_required_value_null_exits(tmp_path: Path) -> None:
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
    assert check_file_extension(allowed, filename) is expected


def test_samplesheet_handle_presets_without_disable_column_uses_cli_flag(monkeypatch: pytest.MonkeyPatch) -> None:
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
    "input_exists,reference,reference_exists,reference_ext,primers,primers_exists,primers_ext,features,features_exists,features_ext,scheduler_valid,expected_min_errors",
    [
        (False, "ref.fa", True, True, "p.fa", True, True, "f.gff", True, True, True, 1),
        (True, "NONE", True, True, "p.fa", True, True, "f.gff", True, True, True, 1),
        (True, "ref.txt", True, False, "p.txt", True, False, "f.txt", True, False, True, 3),
        (True, "ref.fa", False, True, "p.fa", False, True, "f.gff", False, True, False, 4),
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
    expected_min_errors: int,
) -> None:
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
    assert len(errors) >= expected_min_errors


def test_check_sample_sheet_missing_required_columns_exits(monkeypatch: pytest.MonkeyPatch) -> None:
    parser_obj = CLIparser.__new__(CLIparser)
    parser_obj.flags = Namespace(presets=True)

    monkeypatch.setattr("ViroConstrictor.parser.open_sample_sheet", lambda _f: pd.DataFrame({"SAMPLE": ["s1"]}))
    monkeypatch.setattr("ViroConstrictor.parser.check_samplesheet_columns", lambda _df: False)

    with pytest.raises(SystemExit):
        parser_obj._check_sample_sheet("dummy.csv")


def test_check_sample_sheet_empty_returns_empty_df(monkeypatch: pytest.MonkeyPatch) -> None:
    parser_obj = CLIparser.__new__(CLIparser)
    parser_obj.flags = Namespace(presets=True)
    monkeypatch.setattr("ViroConstrictor.parser.open_sample_sheet", lambda _f: pd.DataFrame())

    result = parser_obj._check_sample_sheet("dummy.csv")

    assert result.empty


def test_check_sample_properties_missing_reference_exits() -> None:
    parser_obj = CLIparser.__new__(CLIparser)
    sampleinfo = {"s1": {"REFERENCE": "/tmp/does_not_exist.fasta"}}

    with pytest.raises(SystemExit):
        parser_obj._check_sample_properties(sampleinfo)


def test_check_sample_properties_checks_unique_references(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
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
