from argparse import ArgumentParser
from pathlib import Path

import pandas as pd

from .base_script_class import BaseScript


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

    def __init__(self, input: list[Path | str], output: Path | str) -> None:
        super().__init__(input, output)

    @classmethod
    def add_arguments(cls, parser: ArgumentParser) -> None:
        super().add_arguments(parser)

    def run(self) -> None:
        self._concat_amplicon_covs()

    def _concat_amplicon_covs(self) -> None:
        """
        Concatenates multiple amplicon coverage CSV files into a single CSV file.
        """
        assert isinstance(self.input, list), "Input should be a list of file paths."
        assert isinstance(self.output, (Path, str)), "Output should be a string path for the output file."

        # Read and concatenate input files
        frames = [self._read_csv(file) for file in self.input]
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
