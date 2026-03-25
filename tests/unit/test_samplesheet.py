"""Unit tests for sample sheet discovery helpers."""

from pathlib import Path
from typing import cast

from ViroConstrictor.samplesheet import GetSamples, illumina_sheet, iontorrent_sheet, nanopore_sheet


def test_illumina_sheet_collects_r1_r2_and_ignores_non_fastq(tmp_path: Path) -> None:
    """Test that illumina_sheet correctly collects paired-end reads and ignores non-FASTQ files.

    Verifies that both R1/R2 suffixes are recognized and non-FASTQ files (e.g., .txt) are excluded.

    Parameters
    ----------
    tmp_path : Path
        Temporary directory provided by pytest.
    """
    (tmp_path / "sampleA_R1.fastq").write_text("@r\nAC\n+\n!!\n", encoding="utf-8")
    (tmp_path / "sampleA_R2.fastq.gz").write_text("gz-placeholder", encoding="utf-8")
    (tmp_path / "notes.txt").write_text("ignore", encoding="utf-8")

    result = illumina_sheet(tmp_path)

    assert "sampleA" in result
    assert result["sampleA"]["R1"] == str((tmp_path / "sampleA_R1.fastq").resolve())
    assert result["sampleA"]["R2"] == str((tmp_path / "sampleA_R2.fastq.gz").resolve())


def test_illumina_sheet_supports_dot_separator_and_optional_r(tmp_path: Path) -> None:
    """Test that illumina_sheet supports both underscore and dot separators with optional 'R' prefix.

    Verifies that patterns like sampleB.1.fastq (without explicit R) and sampleB.R2.fastq are both
    correctly recognized and assigned to R1 and R2 respectively.

    Parameters
    ----------
    tmp_path : Path
        Temporary directory provided by pytest.
    """
    (tmp_path / "sampleB.1.fastq").write_text("x", encoding="utf-8")
    (tmp_path / "sampleB.R2.fastq").write_text("x", encoding="utf-8")

    result = illumina_sheet(tmp_path)

    assert result["sampleB"]["R1"] == str((tmp_path / "sampleB.1.fastq").resolve())
    assert result["sampleB"]["R2"] == str((tmp_path / "sampleB.R2.fastq").resolve())


def test_illumina_sheet_nested_directory_and_duplicate_read_overwrites(tmp_path: Path) -> None:
    """Test that illumina_sheet recursively searches nested directories and later matches overwrite earlier ones.

    Verifies that when duplicate sample names with the same read number appear in different directories,
    the later discovered file path overwrites the earlier one.

    Parameters
    ----------
    tmp_path : Path
        Temporary directory provided by pytest.
    """
    nested = tmp_path / "nested"
    nested.mkdir()
    (tmp_path / "sampleC_R1.fastq").write_text("first", encoding="utf-8")
    (nested / "sampleC_R1.fastq").write_text("second", encoding="utf-8")

    result = illumina_sheet(tmp_path)

    assert result["sampleC"]["R1"] == str((nested / "sampleC_R1.fastq").resolve())


def test_nanopore_sheet_collects_first_match_only(tmp_path: Path) -> None:
    """Test that nanopore_sheet returns only the first matched file for each sample name.

    Verifies that when multiple files matching the same sample name are found (e.g., .fastq and .fq),
    only one is retained in the dictionary. Non-FASTQ files like .bam are ignored.

    Parameters
    ----------
    tmp_path : Path
        Temporary directory provided by pytest.
    """
    (tmp_path / "np1.fastq").write_text("x", encoding="utf-8")
    (tmp_path / "np1.fq").write_text("y", encoding="utf-8")
    (tmp_path / "np2.fastq.gz").write_text("z", encoding="utf-8")
    (tmp_path / "not_reads.bam").write_text("ignore", encoding="utf-8")

    result = nanopore_sheet(tmp_path)

    assert result["np1"] in {
        str((tmp_path / "np1.fastq").resolve()),
        str((tmp_path / "np1.fq").resolve()),
    }
    assert result["np2"] == str((tmp_path / "np2.fastq.gz").resolve())
    assert "not_reads" not in result


def test_nanopore_sheet_reads_nested_directories(tmp_path: Path) -> None:
    """Test that nanopore_sheet recursively discovers FASTQ files in nested directories.

    Verifies that nested directory structures are searched and files are resolved to absolute paths.

    Parameters
    ----------
    tmp_path : Path
        Temporary directory provided by pytest.
    """
    nested = tmp_path / "run1"
    nested.mkdir()
    (nested / "np_nested.fastq").write_text("x", encoding="utf-8")

    result = nanopore_sheet(tmp_path)

    assert result == {"np_nested": str((nested / "np_nested.fastq").resolve())}


def test_iontorrent_sheet_matches_expected_extensions(tmp_path: Path) -> None:
    """Test that iontorrent_sheet correctly identifies FASTQ files and ignores non-matching extensions.

    Verifies that both .fastq and .fq extensions (with optional gzip compression) are recognized,
    while unrelated file types like .fast5 are excluded.

    Parameters
    ----------
    tmp_path : Path
        Temporary directory provided by pytest.
    """
    (tmp_path / "ion1.fastq").write_text("x", encoding="utf-8")
    (tmp_path / "ion2.fq.gz").write_text("x", encoding="utf-8")
    (tmp_path / "bad.fast5").write_text("ignore", encoding="utf-8")

    result = iontorrent_sheet(tmp_path)

    assert result["ion1"] == str((tmp_path / "ion1.fastq").resolve())
    assert result["ion2"] == str((tmp_path / "ion2.fq.gz").resolve())
    assert "bad" not in result


def test_getsamples_dispatches_illumina(tmp_path: Path) -> None:
    """Test that GetSamples correctly dispatches to illumina_sheet for the 'illumina' platform.

    Verifies that when platform='illumina' is specified, GetSamples returns a dict with
    R1/R2 keys for paired-end reads.

    Parameters
    ----------
    tmp_path : Path
        Temporary directory provided by pytest.
    """
    (tmp_path / "sampleD_R1.fastq").write_text("x", encoding="utf-8")

    result = GetSamples(tmp_path, "illumina")

    assert isinstance(result, dict)
    sample = cast(dict, result["sampleD"])
    assert sample["R1"] == str((tmp_path / "sampleD_R1.fastq").resolve())


def test_getsamples_dispatches_iontorrent(tmp_path: Path) -> None:
    """Test that GetSamples correctly dispatches to iontorrent_sheet for the 'iontorrent' platform.

    Verifies that when platform='iontorrent' is specified, GetSamples returns a simple
    dict mapping sample names to FASTQ file paths.

    Parameters
    ----------
    tmp_path : Path
        Temporary directory provided by pytest.
    """
    (tmp_path / "sampleE.fastq").write_text("x", encoding="utf-8")

    result = GetSamples(tmp_path, "iontorrent")

    assert result == {"sampleE": str((tmp_path / "sampleE.fastq").resolve())}


def test_getsamples_dispatches_nanopore(tmp_path: Path) -> None:
    """Test that GetSamples correctly dispatches to nanopore_sheet for the 'nanopore' platform.

    Verifies that when platform='nanopore' is specified, GetSamples returns a simple
    dict mapping sample names to FASTQ file paths.

    Parameters
    ----------
    tmp_path : Path
        Temporary directory provided by pytest.
    """
    (tmp_path / "sampleF.fastq.gz").write_text("x", encoding="utf-8")

    result = GetSamples(tmp_path, "nanopore")

    assert result == {"sampleF": str((tmp_path / "sampleF.fastq.gz").resolve())}


def test_getsamples_unknown_platform_returns_empty_dict(tmp_path: Path) -> None:
    """Test that GetSamples returns an empty dictionary for unknown platforms.

    Verifies that when an unrecognized platform string is provided, GetSamples
    returns an empty dict instead of raising an error.

    Parameters
    ----------
    tmp_path : Path
        Temporary directory provided by pytest.
    """
    (tmp_path / "sampleG.fastq").write_text("x", encoding="utf-8")

    result = GetSamples(tmp_path, "unknown")

    assert result == {}
