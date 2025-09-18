from argparse import ArgumentParser
from enum import Enum
from pathlib import Path

import pandas as pd

from helpers.base_script_class import BaseScript  # type: ignore[import]  # noqa: F401,E402


class AltName(Enum):
    ALT = ["alt", "alternative"]

    @staticmethod
    def is_valid_alt_name(string: str) -> bool:
        """
        Checks if the given string is a valid alternative name.
        This should return True in cases like alt1, altleft, ALT, etc.
        """
        return any(string.lower().startswith(alt) for alt in AltName.ALT.value)


class ReadDirection(Enum):
    FORWARD = ["FW", "F", "1", "LEFT", "POSITIVE", "FORWARD", "PLUS"]
    REVERSE = ["RV", "R", "2", "RIGHT", "NEGATIVE", "REVERSE", "MINUS"]

    @staticmethod
    def is_valid_direction(string: str) -> bool:
        """
        Checks if the given string is a valid read direction.
        """
        return any(string.upper() in direction.value for direction in ReadDirection)


class AmpliconCovs(BaseScript):
    """
    Calculates amplicon coverage based on primer and coverage data.

    Parameters
    ----------
    input : str
        Path to the input BED file with primers.
    coverages : str
        Path to the input TSV file with coverages.
    key : str
        Sample ID.
    output : str
        Path to the output CSV file.

    Methods
    -------
    run()
        Executes the amplicon coverage calculation.
    """

    def __init__(
        self, input: Path | str, coverages: Path | str, key: str, output: Path | str
    ) -> None:
        super().__init__(input, output)
        self.coverages = coverages
        self.key = key

    @classmethod
    def add_arguments(cls, parser: ArgumentParser) -> None:
        super().add_arguments(parser)
        parser.add_argument(
            "--coverages",
            metavar="File",
            type=str,
            help="Input TSV file with coverages.",
            required=True,
        )
        parser.add_argument(
            "--key",
            metavar="String",
            type=str,
            help="Sample ID.",
            required=True,
        )

    def run(self) -> None:
        self._calculate_amplicon_coverage()

    def _calculate_amplicon_coverage(self) -> None:
        """
        Calculates amplicon coverage and writes the results to the output file.
        """
        primers = self._open_tsv_file(self.input)
        primers = self._split_primer_names(primers)
        primers["direction"] = primers["direction"].apply(self._standardize_direction)
        amplicon_sizes = self._calculate_amplicon_start_end(primers)

        coverages = self._open_tsv_file(self.coverages, index_col=0)
        amplicon_sizes["coverage"] = amplicon_sizes.apply(
            lambda x: self._calculate_mean_coverage(x, coverages), axis=1
        )
        amplicon_sizes["amplicon_names"] = self._create_amplicon_names_list(primers)

        final_df = pd.DataFrame(
            [amplicon_sizes["coverage"].values],
            columns=amplicon_sizes["amplicon_names"],
            index=[self.key],
        )

        self._write_output(final_df, self.output)

    @staticmethod
    def _open_tsv_file(
        filename: Path | str, index_col: int | None = None
    ) -> pd.DataFrame:
        """
        Opens a TSV file and returns its contents as a pandas DataFrame.
        """
        df = pd.read_csv(filename, sep="\t", header=None, index_col=index_col)
        if df.isnull().to_numpy().any():
            raise ValueError(f"File {filename} contains NaN values.")
        return df

    @staticmethod
    def _split_primer_names(df: pd.DataFrame) -> pd.DataFrame:
        """
        Splits the primer names in the DataFrame into separate columns.
        """

        def _process_primer_row(row: pd.Series) -> pd.Series:
            split_names = row["split_name_list"]

            # we need to enforce the correct dtypes and this is the easiest way
            row["name"] = ""
            row["count"] = 0
            row["alt"] = ""
            row["direction"] = ""

            if len(split_names) == 4:
                if ReadDirection.is_valid_direction(
                    split_names[3]
                ) and AltName.is_valid_alt_name(split_names[2]):
                    row["name"], row["count"], row["alt"], row["direction"] = (
                        split_names
                    )
                elif ReadDirection.is_valid_direction(
                    split_names[2]
                ) and AltName.is_valid_alt_name(split_names[3]):
                    row["name"], row["count"], row["direction"], row["alt"] = (
                        split_names
                    )
                else:
                    raise ValueError(
                        f"Primer name {row[3]} does not match expected format with alt and direction."
                    )
            elif len(split_names) == 3:
                row["name"], row["count"], row["alt"], row["direction"] = (
                    split_names[0],
                    split_names[1],
                    "",
                    split_names[2],
                )
            else:
                raise ValueError(
                    f"Primer name {row[3]} does not contain the expected number of underscores."
                )

            row["count"] = int(row["count"])  # Ensure count is an integer
            row["name"] = str(row["name"])
            row["alt"] = str(row["alt"]) if row["alt"] else ""
            row["direction"] = str(row["direction"])
            return row

        df["split_name_list"] = df[3].str.split("_", expand=False)
        df = df.apply(_process_primer_row, axis=1)
        return df.drop(columns=["split_name_list"])

    @staticmethod
    def _standardize_direction(direction: str) -> str:
        """
        Standardizes the given read direction string to a standardized direction name.
        """

        for read_direction in ReadDirection:
            if direction.upper() in read_direction.value:
                return read_direction.name
        raise ValueError(f"Unrecognized read direction: {direction}")

    @staticmethod
    def _calculate_amplicon_start_end(primers: pd.DataFrame) -> pd.DataFrame:
        """
        Calculates the start and end positions of amplicons based on primer data.
        """
        df = pd.DataFrame(primers["count"].unique(), columns=["amplicon_number"])
        for amplicon_number in df["amplicon_number"]:
            amplicon_group = primers[primers["count"] == amplicon_number]
            df.loc[df["amplicon_number"] == amplicon_number, "start"] = amplicon_group[
                1
            ].min()
            df.loc[df["amplicon_number"] == amplicon_number, "end"] = amplicon_group[
                2
            ].max()
        return df

    @staticmethod
    def _calculate_mean_coverage(
        input_array: pd.Series, coverages: pd.DataFrame
    ) -> float:
        """
        Calculates the mean coverage for a given amplicon.
        """
        return (
            coverages.iloc[int(input_array["start"]) - 1 : int(input_array["end"])]
            .mean()
            .values[0]
        )

    @staticmethod
    def _create_amplicon_names_list(primers: pd.DataFrame) -> list[str]:
        """
        Creates a list of unique amplicon names based on the primers DataFrame.
        """
        amplicon_name = primers.loc[0, "name"]
        return [f"{amplicon_name}_{x}" for x in primers["count"].unique()]

    @staticmethod
    def _write_output(df: pd.DataFrame, output_file: Path | str) -> None:
        """
        Writes the calculated amplicon sizes to a CSV file.
        """
        df.to_csv(output_file, sep=",", index=True, index_label="amplicon_names")


if __name__ == "__main__":
    AmpliconCovs.main()
