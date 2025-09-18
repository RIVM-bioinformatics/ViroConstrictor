"""
Module for handling GenBank files.

This module provides utilities for working with GenBank files, including
checking file extensions, parsing GenBank records, extracting target organisms,
and splitting GenBank files into FASTA and GFF formats.
"""

from pathlib import Path

from BCBio import GFF  # type: ignore
from Bio import SeqIO  # type: ignore


class GenBank:
    """
    A utility class for working with GenBank files.

    This class provides static methods to:
    - Check if a file is a GenBank file based on its extension.
    - Parse GenBank files and extract records.
    - Extract the target organism from GenBank records.
    - Split GenBank files into FASTA and GFF formats.

    Attributes
    ----------
    EXTENSIONS : set[str]
        A set of valid GenBank file extensions.
    """

    EXTENSIONS = {".gb", ".gbk", ".genbank"}

    @staticmethod
    def is_genbank(file_path: Path) -> bool:
        """Check if the file path has a GenBank extension."""
        return any(file_path.suffix == ext for ext in GenBank.EXTENSIONS)

    @staticmethod
    def open_genbank(file_path: Path) -> list[SeqIO.SeqRecord]:
        """Open a GenBank file and return its records."""
        if not GenBank.is_genbank(file_path):
            raise ValueError(f"File {file_path} is not a GenBank file.")
        try:
            return list(SeqIO.parse(file_path, "genbank"))
        except Exception as e:
            raise ValueError(f"Error opening GenBank file: {e}") from e

    @staticmethod
    def _parse_target(records: list[SeqIO.SeqRecord]) -> str:
        """Parse the target organism from GenBank records."""
        organisms: list[str] = [record.annotations.get("organism", "") for record in records]
        organisms = [org.split("(", 1)[0].strip().replace(" ", "_") for org in organisms if org]
        if not all(org == organisms[0] for org in organisms):
            raise ValueError(
                "Not all GenBank records have the same organism annotation.\n"
                "Either edit the GenBank file to have the same organism annotation for all records (strains don't count),\n"
                "or manually provide a target organism name using the --target option."
            )
        return organisms[0]

    @staticmethod
    def split_genbank(file_path: Path) -> tuple[Path, Path, str]:
        """Splits a GenBank file into a reference fasta, a features file and possibly a target file."""

        records = GenBank.open_genbank(file_path)
        with open(file_path.with_suffix(".fasta"), "w", encoding="utf-8") as fasta_file:
            # write all records as fasta to a single file
            for record in records:
                SeqIO.write(record, fasta_file, "fasta")
        fasta_path = file_path.with_suffix(".fasta")

        with open(file_path.with_suffix(".gff"), "w", encoding="utf-8") as gff_file:
            GFF.write(records, gff_file)
        gff_path = file_path.with_suffix(".gff")

        target = GenBank._parse_target(records)

        return fasta_path, gff_path, target
