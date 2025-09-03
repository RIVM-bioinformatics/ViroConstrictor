from argparse import ArgumentParser
from pathlib import Path

import pandas as pd
from Bio import SeqIO

from ..base_script_class import BaseScript


class FilterBestMatchingRef(BaseScript):
    """
    Filters the best matching reference from a BAM file and updates the reference FASTA file.

    Parameters
    ----------
    inputcount : str
        Path to the input CSV file containing mapped reads and reference information.
    inputref : str
        Path to the input FASTA file containing reference sequences.
    filtref : str
        Path to the output FASTA file containing the filtered reference sequence.
    filtcount : str
        Path to the output CSV file containing the filtered reference information.

    Methods
    -------
    run()
        Executes the filtering of the best matching reference.
    """

    def __init__(
        self,
        input: Path | str,
        inputref: Path | str,
        filtref: Path | str,
        output: Path | str,
    ) -> None:
        super().__init__(input, output)
        self.inputref = inputref
        self.filtref = filtref

    @classmethod
    def add_arguments(cls, parser: ArgumentParser) -> None:
        super().add_arguments(parser)
        parser.add_argument(
            "--inputref",
            metavar="File",
            type=str,
            help="Path to the input FASTA file containing reference sequences.",
            required=True,
        )
        parser.add_argument(
            "--filtref",
            metavar="File",
            type=str,
            help="Path to the output FASTA file containing the filtered reference sequence.",
            required=True,
        )

    def run(self) -> None:
        self._filter_best_matching_ref()

    def _filter_best_matching_ref(self) -> None:
        """
        Filters the best matching reference from the input CSV file and updates the reference FASTA file.
        """
        assert isinstance(
            self.input, (Path, str)
        ), "Inputcount should be a string path to the CSV file."
        assert isinstance(
            self.inputref, (Path, str)
        ), "Inputref should be a string path to the FASTA file."
        assert isinstance(
            self.filtref, (Path, str)
        ), "Filtref should be a string path for the output FASTA file."
        assert isinstance(
            self.output, (Path, str)
        ), "Output should be a string path for the output CSV file."

        # Read the input CSV file
        df = pd.read_csv(self.input, index_col=0)

        # Sort by "Mapped Reads" column and keep only the first row
        df = df.sort_values(by="Mapped Reads", ascending=False).iloc[[0]]

        # Get the reference name
        refname = df["Reference"].values[0]

        # Read the reference FASTA file and filter for the reference name
        for record in SeqIO.parse(self.inputref, "fasta"):
            if record.id == refname:
                SeqIO.write(record, self.filtref, "fasta")

        # Write the filtered reference information to the output CSV file
        df.to_csv(self.output)


if __name__ == "__main__":
    FilterBestMatchingRef.main()
