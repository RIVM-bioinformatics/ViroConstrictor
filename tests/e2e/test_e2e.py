from configparser import ConfigParser
from pathlib import Path
from typing import Generator

import pytest

from tests.utils.ena_downloader import download_partial_ena_file
from ViroConstrictor.__main__ import main

LOCK_PATH = Path("tests/e2e/data/output/.snakemake/locks")


def remove_locks(lock_dir: Path) -> None:
    if lock_dir.exists() and lock_dir.is_dir():
        for lock_file in lock_dir.glob("*.lock"):
            lock_file.unlink()


@pytest.fixture(scope="module")
def prepare_files() -> Generator[dict[str, Path], None, None]:
    remove_locks(LOCK_PATH)
    paths = {}

    data_dir = Path(__file__).parent / "data"
    data_dir.mkdir(parents=True, exist_ok=True)

    # input
    filename = "ESIB_EQA_2024_SARS1_01.fastq.gz"
    fastqs_dir = data_dir / "fastqs"
    input_file = fastqs_dir / filename
    if not input_file.exists():
        download_partial_ena_file(
            "SRR9678314",
            fastqs_dir / filename,
            num_reads=1000,
            overwrite=True,
        )
    paths["input"] = fastqs_dir

    # output
    output_dir = data_dir / "output"
    output_dir.mkdir(parents=True, exist_ok=True)
    paths["output"] = output_dir

    # miscellaneous
    paths["reference"] = data_dir / "reference_genome.fasta"
    paths["reference_gb"] = data_dir / "reference_genome.gb"
    paths["reference_MR"] = data_dir / "reference_MR.fasta"
    paths["features"] = data_dir / "ESIB_EQA_2024_SARS1_01_features.gff"
    paths["primers"] = data_dir / "Primers_articv4.1.bed"
    paths["settings"] = data_dir / "test_settings.ini"

    try:
        yield paths
    finally:
        # shutil.rmtree(data_dir)
        pass


def test_main(prepare_files: dict[str, Path]) -> None:
    args = [
        "--input",
        prepare_files["input"].as_posix(),
        "--output",
        prepare_files["output"].as_posix(),
        "--reference",
        prepare_files["reference"].as_posix(),
        "--features",
        prepare_files["features"].as_posix(),
        "--primers",
        prepare_files["primers"].as_posix(),
        "--amplicon-type",
        "fragmented",
        "--platform",
        "nanopore",
        "--target",
        "sars-cov-2",
    ]

    # cwd = os.getcwd()
    # os.rtests/e2e/data/output/.snakemake/locks
    with pytest.raises(SystemExit) as e:
        main(args, settings=prepare_files["settings"].as_posix())
    assert e.value.code == 0, "Main function did not complete successfully"


def test_match_ref(prepare_files: dict[str, Path]) -> None:
    args = [
        "--input",
        prepare_files["input"].as_posix(),
        "--output",
        prepare_files["output"].as_posix(),
        "--reference",
        prepare_files["reference"].as_posix(),
        "--features",
        prepare_files["features"].as_posix(),
        "--primers",
        prepare_files["primers"].as_posix(),
        "--amplicon-type",
        "fragmented",
        "--platform",
        "nanopore",
        "--target",
        "sars-cov-2",
        "-mr",
    ]

    with pytest.raises(SystemExit) as e:
        main(args, settings=prepare_files["settings"].as_posix())
    assert e.value.code == 0, "Main function did not complete successfully"


def test_main_container(prepare_files: dict[str, Path]) -> None:
    args = [
        "--input",
        prepare_files["input"].as_posix(),
        "--output",
        prepare_files["output"].as_posix(),
        "--reference",
        prepare_files["reference"].as_posix(),
        "--features",
        prepare_files["features"].as_posix(),
        "--primers",
        prepare_files["primers"].as_posix(),
        "--amplicon-type",
        "fragmented",
        "--platform",
        "nanopore",
        "--target",
        "sars-cov-2",
    ]

    # Use only repository-local containers to avoid untrusted path input in tests.
    repo_container_folder = (Path(__file__).resolve().parents[2] / "containers").resolve()
    if not any(repo_container_folder.glob("viroconstrictor_*.sif")):
        pytest.skip("No repository-local containers found for container e2e test")
    container_folder = repo_container_folder.as_posix()

    # Update only the required keys via parser instead of raw string rewriting.
    settings_path = prepare_files["settings"]
    original_settings = settings_path.read_text()
    config_reader = ConfigParser()
    config_reader.read_string(original_settings)
    if not config_reader.has_section("REPRODUCTION"):
        config_reader.add_section("REPRODUCTION")
    config_reader["REPRODUCTION"]["repro_method"] = "containers"
    config_reader["REPRODUCTION"]["container_cache_path"] = container_folder
    with settings_path.open("w") as settings_file:
        config_reader.write(settings_file)

    try:
        with pytest.raises(SystemExit) as e:
            main(args, settings=prepare_files["settings"].as_posix())
        assert e.value.code == 0, "Main function did not complete successfully"
    finally:
        # Revert the settings file to its original state
        settings_path.write_text(original_settings)
