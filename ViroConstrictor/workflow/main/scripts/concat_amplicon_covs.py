import os
from argparse import ArgumentParser
from pathlib import Path

import pandas as pd
from helpers.base_script_class import BaseScript  # type: ignore[import]  # noqa: F401,E402


class ConcatAmpliconCovs(BaseScript):
    """
    Concatenates multiple amplicon coverage CSV files into a single CSV file.

    Parameters
    ----------
    input : list[str]
        List of input CSV file paths.
    output : str
        Path to the output CSV file.

    Methods
    -------
    run()
        Executes the concatenation of amplicon coverage files.
    """

    def __init__(self, input: list[Path | str], input_coverages: list[Path | str], output: Path | str) -> None:
        super().__init__(input_coverages, output)
        self.input_coverages = input_coverages

    @classmethod
    def add_arguments(cls, parser: ArgumentParser) -> None:
        super().add_arguments(parser)
        parser.add_argument(
            "--input_coverages",
            metavar="File",
            nargs="+",
            type=str,
            help="List of input CSV file paths.",
            required=True,
        )

    def run(self) -> None:
        self._concat_amplicon_covs()

    def _concat_amplicon_covs(self) -> None:
        """
        Concatenates multiple amplicon coverage CSV files into a single CSV file.
        """
        assert isinstance(self.input, (list, str)), "Input should be a list of file paths."
        assert isinstance(self.output, (Path, str)), "Output should be a string path for the output file."

        if isinstance(self.input, str):
            frames = [self._read_csv(self.input)]
        else:
            frames = [self._read_csv(file) for file in self.input]

        # Read and concatenate input files
        concatenated_frame = pd.concat(frames, sort=True).fillna("NA")

        # Write the concatenated DataFrame to the output file
        concatenated_frame.to_csv(self.output, sep=",", index=True)

    @staticmethod
    def _read_csv(file: Path | str) -> pd.DataFrame:
        """
        Reads a CSV file into a pandas DataFrame.

        Parameters
        ----------
        file : Path | str
            Path to the input CSV file.

        Returns
        -------
        pd.DataFrame
            DataFrame containing the contents of the CSV file.
        """
        df = pd.read_csv(file, sep=",")
        df.set_index("amplicon_names", inplace=True)
        return df


if __name__ == "__main__":
    ConcatAmpliconCovs.main()
