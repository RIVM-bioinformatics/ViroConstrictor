"""
Module for handling GenBank files.

This module provides utilities for working with GenBank files, including
checking file extensions, parsing GenBank records, extracting target organisms,
and splitting GenBank files into FASTA and GFF formats.
"""

import difflib
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
        """Check if the file path has a GenBank extension.
        
        Parameters
        ----------
        file_path : Path
            Path object to check for GenBank file extension.
        
        Returns
        -------
        bool
            True if the file has a valid GenBank extension (.gb, .gbk, .genbank),
            False otherwise.
        """
        return any(file_path.suffix == ext for ext in GenBank.EXTENSIONS)

    @staticmethod
    def open_genbank(file_path: Path) -> list[SeqIO.SeqRecord]:
        """Open a GenBank file and return its records.
        
        Parameters
        ----------
        file_path : Path
            Path to the GenBank file to open.
        
        Returns
        -------
        list[SeqIO.SeqRecord]
            List of parsed SeqIO sequence records from the GenBank file.
        
        Raises
        ------
        ValueError
            If the file is not a GenBank file or if parsing fails.
        """
        if not GenBank.is_genbank(file_path):
            raise ValueError(f"File {file_path} is not a GenBank file.")
        try:
            return list(SeqIO.parse(file_path, "genbank"))
        except Exception as e:
            raise ValueError(f"Error opening GenBank file: {e}") from e

    @staticmethod
    def _parse_target(records: list[SeqIO.SeqRecord]) -> str:
        """Parse the target organism name from GenBank records.
        
        Extracts organism annotations from all records and ensures they are
        sufficiently similar (≥85% sequence similarity) before returning a
        unified target name.
        
        Parameters
        ----------
        records : list[SeqIO.SeqRecord]
            List of SeqIO sequence records from a GenBank file.
        
        Returns
        -------
        str
            The target organism name extracted from record annotations,
            with spaces converted to underscores.
        
        Raises
        ------
        ValueError
            If records have dissimilar organism annotations (similarity < 85%).
        """
        organisms: list[str] = [record.annotations.get("organism", "") for record in records]
        organisms = [org.split("(", 1)[0].strip().replace(" ", "_") for org in organisms if org]

        threshold = 0.85  # similarity threshold
        ref_org = organisms[0]
        if any(difflib.SequenceMatcher(None, ref_org, org).ratio() < threshold for org in organisms):
            raise ValueError(
                "Not all GenBank records have sufficiently similar organism annotations.\n"
                "Either edit the GenBank file to have similar organism annotation for all records (strains don't count),\n"
                "or manually provide a target organism name using the --target option."
            )
        return organisms[0]

    @staticmethod
    def split_genbank(file_path: Path, emit_target: bool = False, output_directory: Path | None = None) -> tuple[Path, Path, str]:
        """Splits a GenBank file into a reference FASTA, a features GFF file, and optionally a target.
        
        Parses the GenBank file and writes its sequences to FASTA format and genomic
        features to GFF3 format in the specified output directory. Optionally extracts
        the target organism name from the GenBank records.
        
        Parameters
        ----------
        file_path : Path
            Path to the GenBank file to split.
        emit_target : bool, optional
            If True, extract and return the target organism name from the GenBank records.
            Default is False.
        output_directory : Path | None, optional
            Directory where split FASTA and GFF files should be written.
            If None, splits are written to the same directory as the source GenBank file.
            Default is None (backward compatibility).
        
        Returns
        -------
        tuple[Path, Path, str]
            A tuple containing:
            - Path to the generated FASTA file
            - Path to the generated GFF file
            - Target organism name (empty string if emit_target is False)
        """

        records = GenBank.open_genbank(file_path)
        output_dir = output_directory if output_directory is not None else file_path.parent
        output_dir.mkdir(parents=True, exist_ok=True)
        fasta_path = output_dir / f"{file_path.stem}.fasta"
        with open(fasta_path, "w", encoding="utf-8") as fasta_file:
            # write all records as fasta to a single file
            for record in records:
                SeqIO.write(record, fasta_file, "fasta")

        gff_path = output_dir / f"{file_path.stem}.gff"
        with open(gff_path, "w", encoding="utf-8") as gff_file:
            GFF.write(records, gff_file)

        target = ""
        if emit_target:
            target = GenBank._parse_target(records)

        return fasta_path, gff_path, target
