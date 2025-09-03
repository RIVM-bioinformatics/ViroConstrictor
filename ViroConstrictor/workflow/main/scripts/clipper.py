from argparse import ArgumentParser
from pathlib import Path

import pysam

from .base_script_class import BaseScript


class Clipper(BaseScript):
    """
    Filters and processes reads from a BAM file based on various criteria.

    Parameters
    ----------
    input : str
        Path to the input BAM file.
    output : str
        Path to the output file where cleaned FASTQ reads will be saved.
    exclude_spliced : bool, optional
        Whether to exclude spliced reads from the output.
    spliced_length_threshold : int, optional
        Threshold for excluding spliced reads based on spliced length.
    min_aligned_length : float, optional
        Minimum length of the aligned part of the read to include.
    max_aligned_length : float, optional
        Maximum length of the aligned part of the read to include.
    only_include_region : str, optional
        Region to include reads where the aligned section starts and ends within the specified region.
    threads : int, optional
        Number of threads for decompressing/compressing the BAM file.

    Methods
    -------
    run()
        Executes the filtering and processing of reads.
    """

    def __init__(
        self,
        input: Path | str,
        output: Path | str,
        exclude_spliced: bool = False,
        spliced_length_threshold: int = 0,
        min_aligned_length: float = 0,
        max_aligned_length: float = 0,
        only_include_region: str | None = None,
        threads: int = 1,
    ) -> None:
        super().__init__(input, output)
        self.exclude_spliced = exclude_spliced
        self.spliced_length_threshold = spliced_length_threshold
        self.min_aligned_length = min_aligned_length
        self.max_aligned_length = max_aligned_length
        self.only_include_region = only_include_region
        self.threads = threads

    @classmethod
    def add_arguments(cls, parser: ArgumentParser) -> None:
        super().add_arguments(parser)
        parser.add_argument(
            "--exclude-spliced",
            action="store_true",
            default=False,
            help="Exclude spliced reads from the output.",
        )
        parser.add_argument(
            "--spliced-length-threshold",
            metavar="INT",
            type=int,
            default=0,
            help="Exclude spliced reads with a spliced length above this threshold.",
        )
        parser.add_argument(
            "--min-aligned-length",
            metavar="FLOAT",
            type=float,
            default=0,
            help="Minimum length of the aligned part of the read to include.",
        )
        parser.add_argument(
            "--max-aligned-length",
            metavar="FLOAT",
            type=float,
            default=0,
            help="Maximum length of the aligned part of the read to include.",
        )
        parser.add_argument(
            "--only-include-region",
            metavar="STRING",
            type=str,
            default=None,
            help="Only include reads where the aligned section starts and ends within the specified region.",
        )
        parser.add_argument(
            "--threads",
            metavar="Number",
            type=int,
            default=1,
            help="Number of threads for decompressing/compressing the BAM file.",
        )

    def run(self) -> None:
        self._process_reads()

    def _process_reads(self) -> None:
        """
        Processes reads from the BAM file based on the provided criteria and writes them to the output file.
        """
        bamfile = pysam.AlignmentFile(self.input, "rb", threads=self.threads)
        reflength = bamfile.lengths[0]
        minimal_read_length = int(reflength * self.min_aligned_length)
        maximum_read_length = int(
            reflength if self.max_aligned_length == 0 else self.max_aligned_length
        )

        include_region_start = (
            int(self.only_include_region.split(":")[0])
            if self.only_include_region
            else None
        )
        include_region_end = (
            int(self.only_include_region.split(":")[1])
            if self.only_include_region
            else None
        )

        with open(self.output, "w") as fileout:
            for read in bamfile:
                read_start = read.query_alignment_start
                read_end = read.query_alignment_end
                ref_start = read.reference_start
                ref_end = read.reference_end
                trimmed_seq = read.query_alignment_sequence
                trimmed_qual = read.qual[read_start:read_end]

                if read.is_reverse:
                    complement = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
                    bases = [complement[base] for base in trimmed_seq]
                    trimmed_seq = "".join(bases[::-1])
                    trimmed_qual = trimmed_qual[::-1]

                if len(trimmed_seq) == 0:
                    continue

                if self.exclude_spliced:
                    cigartuples = self._split_cigar(read.cigarstring)
                    if self._is_spliced(cigartuples) and (
                        self._get_largest_spliced_len(cigartuples)
                        > self.spliced_length_threshold
                    ):
                        continue

                if read.query_alignment_length <= minimal_read_length:
                    continue

                if read.query_alignment_length >= maximum_read_length:
                    continue

                if include_region_start is not None and include_region_end is not None:
                    if ref_start < include_region_start or ref_end > include_region_end:
                        continue

                fileout.write(f"@{read.query_name}\n{trimmed_seq}\n+\n{trimmed_qual}\n")

    @staticmethod
    def _split_cigar(cigar: str) -> list[tuple[int, str]]:
        """
        Split a CIGAR string into a list of tuples.

        Parameters
        ----------
        cigar : str
            The CIGAR string to be split.

        Returns
        -------
        list
            A list of tuples, where each tuple contains two elements:
            - The number of consecutive characters.
            - The character itself.
        """
        cigar_tuples = []
        current_number = ""
        for char in cigar:
            if char.isdigit():
                current_number += char
            else:
                cigar_tuples.append((int(current_number), char))
                current_number = ""
        return cigar_tuples

    @staticmethod
    def _is_spliced(cigar_tuples: list[tuple[int, str]]) -> bool:
        """
        Check if the given list of cigar tuples contains any spliced regions.

        Parameters
        ----------
        cigar_tuples : list
            A list of cigar tuples representing the alignment.

        Returns
        -------
        bool
            True if there are spliced regions, False otherwise.
        """
        return any(cigar_tuple[1] == "N" for cigar_tuple in cigar_tuples)

    @staticmethod
    def _get_largest_spliced_len(cigar_tuples: list[tuple[int, str]]) -> int:
        """
        Calculate the largest spliced length from a list of cigar tuples.

        Parameters
        ----------
        cigar_tuples : list
            A list of cigar tuples representing the cigar string.

        Returns
        -------
        int
            The largest spliced length found in the cigar tuples.
        """
        largest_spliced_len = 0
        current_spliced_len = 0
        for cigar_tuple in cigar_tuples:
            if cigar_tuple[1] == "N":
                current_spliced_len += cigar_tuple[0]
            else:
                largest_spliced_len = max(largest_spliced_len, current_spliced_len)
                current_spliced_len = 0
        return largest_spliced_len


if __name__ == "__main__":
    Clipper.main()
