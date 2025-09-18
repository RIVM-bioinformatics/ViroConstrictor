from argparse import ArgumentParser
from pathlib import Path

import pandas as pd
import pysam
from helpers.base_script_class import BaseScript  # type: ignore[import]  # noqa: F401,E402


class CountMappedReads(BaseScript):
    """
    Counts mapped reads, mismatches, and sequence identity for each reference in a BAM file.

    Parameters
    ----------
    input : str
        Path to the input BAM file.
    output : str
        Path to the output CSV file.

    Methods
    -------
    run()
        Executes the counting of mapped reads and writes the results to the output file.
    """

    def __init__(self, input: Path | str, output: Path | str) -> None:
        super().__init__(input, output)

    @classmethod
    def add_arguments(cls, parser: ArgumentParser) -> None:
        super().add_arguments(parser)

    def run(self) -> None:
        self._count_mapped_reads()

    def _count_mapped_reads(self) -> None:
        """
        Counts mapped reads, mismatches, and sequence identity for each reference in the BAM file.
        """
        assert isinstance(self.input, (Path, str)), "Input should be a string path to the BAM file."
        assert isinstance(self.output, (Path, str)), "Output should be a string path for the output CSV file."

        # Open BAM file
        bamfile = pysam.AlignmentFile(self.input, "rb")

        # Loop over reference headers and count mapped reads, mismatches, and sequence identity
        data = []
        for ref in bamfile.references:
            mapped_reads = 0
            total_mismatches = 0
            total_identity = 0
            for read in bamfile.fetch(ref):
                if not read.is_unmapped:
                    mapped_reads += 1
                    total_mismatches += read.get_tag("NM")
                    total_identity += 1 - (read.get_tag("NM") / read.query_alignment_length)
            if mapped_reads > 0:
                avg_mismatches = total_mismatches / mapped_reads
                avg_identity = total_identity / mapped_reads
            else:
                avg_mismatches = 0
                avg_identity = 0
            data.append(
                {
                    "Reference": ref,
                    "Mapped Reads": mapped_reads,
                    "Avg. Mismatches per Read": avg_mismatches,
                    "Avg. Sequence Identity": avg_identity,
                }
            )

        # Create pandas DataFrame
        df = pd.DataFrame(data)

        # Write output to CSV file
        df.to_csv(self.output, index=False)


if __name__ == "__main__":
    CountMappedReads.main()
