from configparser import ConfigParser
from pathlib import Path
from typing import Generator

import pytest

from tests.utils.ena_downloader import download_partial_ena_file
from ViroConstrictor.__main__ import main


@pytest.fixture(scope="module")
def prepare_files() -> Generator[dict[str, Path], None, None]:
    paths = {}

    data_dir = Path(__file__).parent / "data"
    data_dir.mkdir(parents=True, exist_ok=True)

    # input
    filename = "test_data_illumina_miseq_h1n1.fastq"
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
    paths["reference"] = data_dir / "test_reference.fasta"
    paths["reference_gb"] = data_dir / "test_reference.gb"
    paths["features"] = data_dir / "test_features.gff"
    paths["primers"] = data_dir / "test_primers.bed"
    paths["settings"] = data_dir / "test_settings.ini"

    try:
        yield paths
    finally:
        # shutil.rmtree(data_dir)
        pass


# def test_main(prepare_files: dict[str, Path]) -> None:
#     args = [
#         "--input",
#         prepare_files["input"].as_posix(),
#         "--output",
#         prepare_files["output"].as_posix(),
#         "--reference",
#         prepare_files["reference"].as_posix(),
#         "--features",
#         prepare_files["features"].as_posix(),
#         "--primers",
#         prepare_files["primers"].as_posix(),
#         "--amplicon-type",
#         "fragmented",
#         "--platform",
#         "nanopore",
#         "--target",
#         "Influenza_A",
#     ]

#     # cwd = os.getcwd()
#     # os.rtests/e2e/data/output/.snakemake/locks
#     with pytest.raises(SystemExit) as e:
#         main(args, settings=prepare_files["settings"].as_posix())
#         assert e.value.code == 0, "Main function did not complete successfully"


# def test_main_container(prepare_files: dict[str, Path]) -> None:
#     args = [
#         "--input",
#         prepare_files["input"].as_posix(),
#         "--output",
#         prepare_files["output"].as_posix(),
#         "--reference",
#         prepare_files["reference"].as_posix(),
#         "--features",
#         prepare_files["features"].as_posix(),
#         "--primers",
#         prepare_files["primers"].as_posix(),
#         "--amplicon-type",
#         "fragmented",
#         "--platform",
#         "nanopore",
#         "--target",
#         "Influenza_A",
#     ]

#     user_path = Path("~/.ViroConstrictor_defaultprofile.ini").expanduser()
#     if user_path.exists():
#         config_reader = ConfigParser()
#         config_reader.read(user_path)
#         container_folder = config_reader["REPRODUCTION"]["container_cache_path"]
#     else:
#         container_folder = Path("./containers").as_posix()

#     # Read and update the settings file dynamically
#     settings_lines = prepare_files["settings"].read_text().splitlines()
#     updated_settings = []
#     for line in settings_lines:
#         if line.startswith("repro_method"):
#             updated_settings.append("repro_method = containers")
#         elif line.startswith("container_cache_path"):
#             updated_settings.append(f"container_cache_path = {container_folder}")
#         else:
#             updated_settings.append(line)
#     prepare_files["settings"].write_text("\n".join(updated_settings))

#     try:
#         with pytest.raises(SystemExit) as e:
#             main(args, settings=prepare_files["settings"].as_posix())
#             assert e.value.code == 0, "Main function did not complete successfully"
#     finally:
#         # Revert the settings file to its original state
#         prepare_files["settings"].write_text("\n".join(settings_lines))


def test_main_genbank(prepare_files: dict[str, Path]) -> None:
    args = [
        "--input",
        prepare_files["input"].as_posix(),
        "--output",
        prepare_files["output"].as_posix(),
        "--reference",
        prepare_files["reference_gb"].as_posix(),
        "--primers",
        prepare_files["primers"].as_posix(),
        "--amplicon-type",
        "fragmented",
        "--platform",
        "nanopore",
    ]

    with pytest.raises(SystemExit) as e:
        main(args, settings=prepare_files["settings"].as_posix())
        assert e.value.code == 0, "Main function did not complete successfully"
