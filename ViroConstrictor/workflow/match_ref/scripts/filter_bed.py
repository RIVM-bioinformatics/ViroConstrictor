import os
from argparse import ArgumentParser
from pathlib import Path
from typing import TYPE_CHECKING

import pandas as pd
from helpers.base_script_class import BaseScript  # type: ignore[import]  # noqa: F401,E402


class FilterBed(BaseScript):
    """
    Filters a BED file based on reference data and updates the statistics.

    Parameters
    ----------
    input : str
        Path to the input BED file.
    refdata : str
        Path to the reference data CSV file.
    output : str
        Path to the filtered BED file.
    updatedstats : str
        Path to the updated statistics CSV file.

    Methods
    -------
    run()
        Executes the filtering of the BED file and updates the statistics.
    """

    def __init__(
        self,
        input: Path | str,
        refdata: Path | str,
        output: Path | str,
        updatedstats: Path | str,
    ) -> None:
        super().__init__(input, output)
        self.refdata = refdata
        self.updatedstats = updatedstats

    @classmethod
    def add_arguments(cls, parser: ArgumentParser) -> None:
        super().add_arguments(parser)
        parser.add_argument(
            "--refdata",
            metavar="File",
            type=str,
            help="Path to the reference data CSV file.",
            required=True,
        )
        parser.add_argument(
            "--updatedstats",
            metavar="File",
            type=str,
            help="Path to the updated statistics CSV file.",
            required=True,
        )

    def run(self) -> None:
        self._filter_bed()

    def _filter_bed(self) -> None:
        """
        Filters the BED file based on reference data and updates the statistics.
        """
        assert isinstance(self.input, (Path, str)), "Input should be a string path to the BED file."
        assert isinstance(self.refdata, (Path, str)), "Reference data should be a string path to the CSV file."
        assert isinstance(self.output, (Path, str)), "Output should be a string path for the filtered BED file."
        assert isinstance(self.updatedstats, (Path, str)), "Updated stats should be a string path for the CSV file."

        # Read reference data
        ref_df = pd.read_csv(self.refdata, keep_default_na=False)

        # Read BED file
        bed_df = pd.read_csv(
            self.input,
            sep="\t",
            names=["ref", "start", "stop", "name", "score", "strand"],
            index_col=False,
        )

        # Filter BED file based on reference data
        bed_df = bed_df.loc[bed_df["ref"].isin(ref_df["seqrecord_name"].tolist())].reset_index(drop=True)

        # Replace "ref" column values in BED file with corresponding "Reference" column values from reference data
        bed_df["ref"] = bed_df["ref"].map(ref_df.set_index("seqrecord_name")["Reference"])

        # Write filtered BED file to output
        bed_df.to_csv(self.output, sep="\t", index=False, header=False)

        # Update statistics and write to updated stats file
        ref_df["Primer_file"] = Path(self.output).resolve()
        ref_df.to_csv(self.updatedstats, index=False)


if __name__ == "__main__":
    FilterBed.main()
