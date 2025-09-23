"""Module for generating FASTA files for testing purposes."""

from pathlib import Path
from random import choice


def generate_fasta_file(
    output_path: Path,
    length: int,
    number_of_sequences: int = 1,
    descriptions: list[str] = ["Test sequence"],
) -> None:
    """
    Generate a FASTA file with a specified number of sequences and length.

    Parameters
    ----------
    output_path : Path
        Path to the output FASTA file.
    length : int
        Length of each sequence.
    number_of_sequences : int, optional
        Number of sequences to generate (default is 1).
    description : str, optional
        Description for each sequence (default is "Test sequence").

    Returns
    -------
    None
        Writes the generated sequences to the specified FASTA file.
    """
    with open(output_path, "w", encoding="utf-8") as fasta_file:
        for seq_num, desc in zip(range(number_of_sequences), descriptions):
            sequence = "".join(choice("ACGT") for _ in range(length))
            fasta_file.write(f">{desc} {seq_num + 1}\n{sequence}\n")
