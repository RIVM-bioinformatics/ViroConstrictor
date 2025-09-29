"""
Script to download a partial ENA file from the European Nucleotide Archive (ENA) using the ENA API.
Useful for testing purposes, since the full file can be large,
and now we don't have to store it in the repo.
"""

import gzip
from pathlib import Path

import requests

ACCESSION = "SRR30230972"
OUTPUT_PATH = "tests/e2e/data/fastqs/test_data_nanopore_h1n1.fastq"
READS = 100


def validate_output_path(output_path: str | Path, overwrite: bool) -> Path:
    """
    Validate the output path and create the directory if it doesn't exist.
    """
    if isinstance(output_path, str):
        output_path = Path(output_path)
    if output_path.suffix != ".fastq":
        raise ValueError("Output file must have .fastq extension")
    if output_path.exists() and not overwrite:
        raise FileExistsError(f"File {output_path} already exists")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    return output_path


def fetch_fastq_url(accession_num: str) -> str:
    """
    The ENA API returns a TSV file with the fastq_ftp field.
    The fastq_ftp field contains a list of URLs separated by semicolons.
    We only need the first URL, so we split the string by semicolons and take the first one.
    """
    url = f"https://www.ebi.ac.uk/ena/portal/api/filereport?accession={accession_num}" "&result=read_run&fields=fastq_ftp&format=tsv&download=true"

    response = requests.get(url, stream=False, timeout=10)
    if response.status_code != 200:
        raise requests.exceptions.HTTPError(f"Failed to download file: {response.status_code} {response.reason}")
    lines = response.text.strip().split("\n")
    if len(lines) < 2:
        raise ValueError("Invalid TSV response format")

    header = lines[0].split("\t")
    values = lines[1].split("\t")

    if "fastq_ftp" not in header:
        raise ValueError("fastq_ftp field not found in response")

    fastq_ftp_index = header.index("fastq_ftp")
    fastq_url = values[fastq_ftp_index].split(";")[0]

    if not fastq_url:
        raise ValueError("No fastq URLs found")
    return fastq_url


def download_fastq_file(fastq_url: str) -> requests.Response:
    """
    Download the fastq file from the ENA FTP server,
    using the URL obtained from the ENA API.
    """
    fastq_url = f"https://{fastq_url}"

    fastq_response = requests.get(fastq_url, stream=True, timeout=10)
    if fastq_response.status_code != 200:
        raise requests.exceptions.HTTPError(f"Failed to download fastq file: {fastq_response.status_code} {fastq_response.reason}")

    return fastq_response


def write_fastq_reads_to_file(num_reads: int, output_file: Path, fastq_response: requests.Response) -> None:
    """
    Write the reads from the fastq response to a file.
    The fastq response is a gzip file, so we need to decompress it.
    We only write the amount or reads equal to num_reads to the file.
    """
    with open(output_file, "wb") as f:
        read_counter = 0
        line_buffer: list[str] = []

        for line in gzip.open(fastq_response.raw, "rt"):
            line_buffer.append(line)
            if len(line_buffer) == 4:
                f.write("".join(line_buffer).encode("utf-8"))
                line_buffer = []
                read_counter += 1
                if read_counter >= num_reads:
                    break


def download_partial_ena_file(
    accession_num: str,
    output_path: str | Path,
    num_reads: int = 1_000,
    overwrite: bool = False,
) -> None:
    """
    Download X amount of reads from a fastq file from the ENA FTP server.
    """
    output_file = validate_output_path(output_path, overwrite)
    fastq_url = fetch_fastq_url(accession_num)
    fastq_response = download_fastq_file(fastq_url)
    write_fastq_reads_to_file(num_reads, output_file, fastq_response)


if __name__ == "__main__":
    download_partial_ena_file(ACCESSION, OUTPUT_PATH, READS, overwrite=True)
