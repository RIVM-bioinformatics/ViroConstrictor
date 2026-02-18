"""
Amplicon Coverage Calculator for ViroConstrictor.

This module calculates amplicon coverage based on primer and coverage data for viral genome
sequencing analysis. It provides flexible primer name parsing to handle various naming
conventions and computes mean coverage across amplicon regions.

The script processes BED files containing primer information and TSV files with coverage
data to produce amplicon-specific coverage statistics. It supports complex primer naming
schemes including alternative primers and various direction indicators.

Classes
-------
AltName : Enum
    Enum to represent alternative primer types.
ReadDirection : Enum
    Enum to represent primer read directions (forward/reverse).
PrimerInfo : dataclass
    Dataclass to hold parsed primer information.
PrimerNameParser : class
    Flexible parser for primer names with regex pattern matching.
AmpliconCovs : class
    Main class that calculates amplicon coverage from primer and coverage data.

Functions
---------
The module can be executed as a script with command-line arguments:
    python amplicon_covs.py --input primers.bed --coverages coverage.tsv --key sample_id --output results.csv

Examples
--------
>>> parser = PrimerNameParser()
>>> primer_info = parser.parse("ncov-2019_1_LEFT")
>>> print(primer_info.name, primer_info.count, primer_info.direction.name)
ncov-2019 1 FORWARD

>>> covs = AmpliconCovs("primers.bed", "coverage.tsv", "sample1", "output.csv")
>>> covs.run()

Notes
-----
Supported primer name formats:
- name_number_direction (e.g., "ncov-2019_1_LEFT")
- name_number_alt_direction (e.g., "ncov-2019_1_alt_LEFT")
- name_number_direction_alt (e.g., "ncov-2019_1_LEFT_alt1")
- Complex alt formats (e.g., "HAV_1_alt-IB_LEFT")
- Multiple underscores in names (e.g., "virus_strain_1_LEFT")
- Multiple underscores with alt (e.g., "virus_strain_1_alt_LEFT")
- name_number-alt_direction (e.g., "MuV-NGS_19-alt22_LEFT

Valid direction indicators:
    FW, F, LEFT, POSITIVE, FORWARD, PLUS (forward)
    RV, R, RIGHT, NEGATIVE, REVERSE, MINUS (reverse)

"""

import re
from argparse import ArgumentParser
from dataclasses import dataclass
from enum import Enum
from pathlib import Path

import pandas as pd
from helpers.base_script_class import BaseScript  # type: ignore[import]  # noqa: F401,E402


class AltName(Enum):
    """Enum to represent alternative primers."""

    ALT = ["alt", "alternative"]

    @staticmethod
    def is_valid_alt_name(string: str) -> bool:
        """
        Checks if the given string is a valid alternative name.
        This should return True in cases like alt1, altleft, ALT, etc.
        """
        return any(string.lower().startswith(alt) for alt in AltName.ALT.value)


class ReadDirection(Enum):
    """Enum to represent read directions."""

    FORWARD = ["FW", "F", "LEFT", "POSITIVE", "FORWARD", "PLUS"]
    REVERSE = ["RV", "R", "RIGHT", "NEGATIVE", "REVERSE", "MINUS"]

    @staticmethod
    def is_valid_direction(string: str) -> bool:
        """
        Checks if the given string is a valid read direction.
        """
        return any(string.upper() in direction.value for direction in ReadDirection)

    @staticmethod
    def from_string(string: str) -> "ReadDirection":
        """
        Converts a string to a ReadDirection enum member.
        Raises ValueError if the string does not correspond to any ReadDirection.
        """
        for direction in ReadDirection:
            if string.upper() in direction.value:
                return direction
        raise ValueError(f"Unrecognized read direction: {string}")


@dataclass
class PrimerInfo:
    """Dataclass to hold parsed primer information."""

    name: str
    count: int
    alt: bool
    direction: ReadDirection
    original_string: str


class PrimerNameParser:
    """Class to parse and validate primer names."""

    def __init__(self):
        # Regex patterns for different primer name formats
        self.patterns = [
            # Pattern 1: name_number_direction (e.g., "ncov-2019_1_LEFT")
            r"^([^_]+)_(\d+)_([^_]+)$",
            # Pattern 2: name_number_alt_direction (e.g., "ncov-2019_1_alt_LEFT")
            r"^([^_]+)_(\d+)_(alt\w*)_([^_]+)$",
            # Pattern 3: name_number_direction_alt (e.g., "ncov-2019_1_LEFT_alt1")
            r"^([^_]+)_(\d+)_([^_]+)_(alt\w*)$",
            # Pattern 4: More complex alt formats (e.g., "HAV_1_alt-IB_LEFT")
            r"^([^_]+)_(\d+)_(alt[^_]*)_([^_]+)$",
            # Pattern 5: Multiple underscores in name (e.g., "virus_strain_1_LEFT")
            r"^(.+)_(\d+)_([^_]+)$",
            # Pattern 6: Multiple underscores with alt (e.g., "virus_strain_1_alt_LEFT")
            r"^(.+)_(\d+)_(alt\w*)_([^_]+)$",
            # Pattern 7: name_number-alt_direction (e.g., "MuV-NGS_19-alt22_LEFT")
            r"^([^_]+)_(\d+)-(alt\w*)_([^_]+)$",
        ]

    def parse(self, primer_name: str) -> PrimerInfo:
        """Parses the primer name and returns a PrimerInfo dataclass."""
        cleaned_name = primer_name.lstrip(">")  # Remove leading '>' if present
        for pattern in self.patterns:
            match = re.match(pattern, cleaned_name, re.IGNORECASE)

            if match:
                return self._extract_info(match, cleaned_name)

        return self._fallback_parse(cleaned_name)

    def _extract_info(self, match: re.Match[str], original_string: str) -> PrimerInfo:
        groups: tuple[str, ...] = match.groups()

        alt = False
        read_direction, count = None, None
        i, j, k = None, None, None

        for i, group in enumerate(groups):
            if AltName.is_valid_alt_name(group):
                alt = True
                break
        if not alt:
            i = -1  # so that it does not interfere with name extraction

        for j, group in enumerate(groups):
            if ReadDirection.is_valid_direction(group):
                read_direction = ReadDirection.from_string(group)
                break
        if read_direction is None:
            raise ValueError(f"Unrecognized read direction in primer name: {original_string}")

        for k, group in enumerate(groups):
            if group.isdigit():
                count = int(group)
                break
        if count is None:
            raise ValueError(f"Primer number not found in primer name: {original_string}")

        name_parts = [group for idx, group in enumerate(groups) if idx not in (i, j, k)]

        name = "_".join(name_parts) if name_parts else "unknown"

        return PrimerInfo(name=name, count=count, alt=alt, direction=read_direction, original_string=original_string)

    def _fallback_parse(self, primer_name: str) -> PrimerInfo:

        name, count, alt, read_direction = None, None, None, None
        parts = primer_name.split("_")
        if len(parts) <= 2:
            raise ValueError(f"Primer name {primer_name} does not match expected formats.")

        if len(parts) == 3:
            name, count_str, direction_str = parts
            return PrimerInfo(
                name=name,
                count=int(count_str),
                alt=None,
                direction=ReadDirection.from_string(direction_str),
                original_string=primer_name,
            )
        elif len(parts) == 4:
            name, count, part3, part4 = parts
            if AltName.is_valid_alt_name(part3):
                alt = AltName.ALT
                read_direction = ReadDirection.from_string(part4)
            elif AltName.is_valid_alt_name(part4):
                alt = AltName.ALT
                read_direction = ReadDirection.from_string(part3)
            else:
                raise ValueError(f"Primer name {primer_name} does not match expected formats.")
            return PrimerInfo(
                name=name,
                count=int(count),
                alt=alt,
                direction=read_direction,
                original_string=primer_name,
            )
        else:
            raise ValueError(f"Primer name {primer_name} does not match expected formats.")


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

    def __init__(self, input: Path | str, coverages: Path | str, key: str, output: Path | str) -> None:
        super().__init__(input, output)
        self.coverages = coverages
        self.key = key
        self.parser = PrimerNameParser()  # Initialize the parser

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
        primers = self._open_tsv_file(self.input)  # _open_tsv_file handles empty files, will return empty df
        if primers.empty:
            final_df = pd.DataFrame(columns=[], index=[self.key])
            self._write_output(final_df, self.output)
            return
        primers = self._split_primer_names(primers)
        amplicon_sizes = self._calculate_amplicon_start_end(primers)

        coverages = self._open_tsv_file(self.coverages, index_col=0)
        amplicon_sizes["coverage"] = amplicon_sizes.apply(lambda x: self._calculate_mean_coverage(x, coverages), axis=1)
        amplicon_sizes["amplicon_names"] = self._create_amplicon_names_list(primers)

        final_df = pd.DataFrame(
            [amplicon_sizes["coverage"].values],
            columns=amplicon_sizes["amplicon_names"],
            index=[self.key],
        )

        self._write_output(final_df, self.output)

    @staticmethod
    def _open_tsv_file(filename: Path | str, index_col: int | None = None) -> pd.DataFrame:
        """
        Opens a TSV file and returns its contents as a pandas DataFrame.
        """
        try:
            df = pd.read_csv(filename, sep="\t", header=None, index_col=index_col, keep_default_na=False)
        except pd.errors.EmptyDataError:
            return pd.DataFrame()

        if df.isnull().to_numpy().any():
            raise ValueError(f"File {filename} contains NaN values.")
        return df

    def _split_primer_names(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Splits the primer names in the DataFrame into separate columns using the parser.
        """

        def _process_primer_row(row: pd.Series) -> pd.Series:
            primer_name = row[3]  # Get the primer name from column 3

            # Use the parser to parse the primer name
            parser_primer_name = self.parser.parse(primer_name)

            # Set the parsed values using attribute access
            row["name"] = parser_primer_name.name
            row["count"] = int(parser_primer_name.count)
            row["alt"] = parser_primer_name.alt
            row["direction"] = parser_primer_name.direction

            return row

        # Apply the parsing function to each row
        df = df.apply(_process_primer_row, axis=1)
        return df

    # TODO: reduce complexity of this function
    @staticmethod
    def _calculate_amplicon_start_end(primers: pd.DataFrame) -> pd.DataFrame:
        """
        Calculates non-overlapping start and end positions of amplicons based on primer data.

        The amplicon regions are partitioned to avoid overlaps:
        - Start: end position of the PREVIOUS amplicon's RIGHT primer
            (or current LEFT primer's end for first amplicon)
        - End: start position of the NEXT amplicon's LEFT primer
            (or current RIGHT primer's start for last amplicon)

        This ensures true non-overlapping amplicon regions for accurate coverage calculation.
        """
        amplicon_numbers = sorted(primers["count"].unique())
        df = pd.DataFrame(amplicon_numbers, columns=["amplicon_number"])

        for idx, amplicon_number in enumerate(amplicon_numbers):
            amplicon_group = primers[primers["count"] == amplicon_number]

            # Get forward (LEFT) and reverse (RIGHT) primers for this amplicon
            forward_primers = amplicon_group[amplicon_group["direction"] == ReadDirection.FORWARD]
            reverse_primers = amplicon_group[amplicon_group["direction"] == ReadDirection.REVERSE]

            # Determine start position
            if idx == 0:
                # First amplicon: start at end of own LEFT primer
                if not forward_primers.empty:
                    amplicon_start = forward_primers[2].max()  # Column 2 is the end position
                else:
                    amplicon_start = amplicon_group[1].min()
            else:
                # Subsequent amplicons: start at end of previous amplicon's RIGHT primer
                prev_amplicon_number = amplicon_numbers[idx - 1]
                prev_amplicon_group = primers[primers["count"] == prev_amplicon_number]
                prev_reverse_primers = prev_amplicon_group[prev_amplicon_group["direction"] == ReadDirection.REVERSE]

                if not prev_reverse_primers.empty:
                    amplicon_start = prev_reverse_primers[2].max()  # Column 2 is the end position
                else:
                    # Fallback to own LEFT primer end
                    amplicon_start = forward_primers[2].max() if not forward_primers.empty else amplicon_group[1].min()

            # Determine end position
            if idx == len(amplicon_numbers) - 1:
                # Last amplicon: end at start of own RIGHT primer
                if not reverse_primers.empty:
                    amplicon_end = reverse_primers[1].min()  # Column 1 is the start position
                else:
                    amplicon_end = amplicon_group[2].max()
            else:
                # Non-last amplicons: end at start of next amplicon's LEFT primer
                next_amplicon_number = amplicon_numbers[idx + 1]
                next_amplicon_group = primers[primers["count"] == next_amplicon_number]
                next_forward_primers = next_amplicon_group[next_amplicon_group["direction"] == ReadDirection.FORWARD]

                if not next_forward_primers.empty:
                    amplicon_end = next_forward_primers[1].min()  # Column 1 is the start position
                else:
                    # Fallback to own RIGHT primer start
                    amplicon_end = reverse_primers[1].min() if not reverse_primers.empty else amplicon_group[2].max()

            # Validate that start < end
            if amplicon_start >= amplicon_end:
                raise ValueError(
                    f"Invalid amplicon {amplicon_number}: start position ({amplicon_start}) "
                    f"is greater than or equal to end position ({amplicon_end}). "
                    f"This may indicate malformed primer data."
                )

            df.loc[df["amplicon_number"] == amplicon_number, "start"] = amplicon_start
            df.loc[df["amplicon_number"] == amplicon_number, "end"] = amplicon_end

        return df

    @staticmethod
    def _calculate_mean_coverage(input_array: pd.Series, coverages: pd.DataFrame) -> float:
        """
        Calculates the mean coverage for a given amplicon.
        """
        return round(float(coverages.iloc[int(input_array["start"]) - 1 : int(input_array["end"])].mean().values[0]), 2)

    @staticmethod
    def _create_amplicon_names_list(primers: pd.DataFrame) -> list[str]:
        """
        Creates a list of unique amplicon names based on the primers DataFrame.
        Amplicon numbers are zero-padded to 3 digits for proper sorting.
        """
        amplicon_name = primers.loc[0, "name"]
        return [f"{amplicon_name}_{str(x).zfill(3)}" for x in sorted(primers["count"].unique())]

    @staticmethod
    def _write_output(df: pd.DataFrame, output_file: Path | str) -> None:
        """
        Writes the calculated amplicon sizes to a CSV file.
        """
        df.to_csv(output_file, sep=",", index=True, index_label="amplicon_names")


if __name__ == "__main__":
    AmpliconCovs.main()
