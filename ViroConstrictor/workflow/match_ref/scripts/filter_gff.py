from argparse import ArgumentParser
from pathlib import Path

import pandas as pd
from AminoExtract.logging import log as AminoExtractLogger
from AminoExtract.reader import SequenceReader
from helpers.base_script_class import BaseScript  # type: ignore[import]  # noqa: F401,E402


class FilterGff(BaseScript):
    """
    Filters a GFF file based on reference data and updates the statistics.

    Parameters
    ----------
    input : str
        Path to the input GFF file.
    refdata : str
        Path to the reference data CSV file.
    output : str
        Path to the filtered GFF file.
    updatedstats : str
        Path to the updated statistics CSV file.

    Methods
    -------
    run()
        Executes the filtering of the GFF file and updates the statistics.
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
        self._filter_gff()

    def _filter_gff(self) -> None:
        """
        Filters the GFF file based on reference data and updates the statistics.
        """
        assert isinstance(self.input, (Path, str)), "Input should be a string path to the GFF file."
        assert isinstance(self.refdata, (Path, str)), "Reference data should be a string path to the CSV file."
        assert isinstance(self.output, (Path, str)), "Output should be a string path for the filtered GFF file."
        assert isinstance(self.updatedstats, (Path, str)), "Updated stats should be a string path for the CSV file."

        # Read the GFF file
        reader = SequenceReader(logger=AminoExtractLogger, verbose=False)
        gff = reader.read_gff(file=self.input)

        # Read the reference data
        df = pd.read_csv(self.refdata, keep_default_na=False)

        # Filter the GFF file based on reference data
        seqrecord_names = df["seqrecord_name"].tolist()
        gff.df = gff.df.loc[gff.df["seqid"].isin(seqrecord_names)].reset_index(drop=True)

        # Update the "seqid" column in the GFF file with the "Reference" column from the reference data
        gff.df["seqid"] = gff.df["seqid"].map(df.set_index("seqrecord_name")["Reference"])

        # Write the filtered GFF file to the output
        gff.export_gff_to_file(self.output)

        # Update the statistics and write to the updated stats file
        df["Feat_file"] = Path(self.output).resolve()
        df.to_csv(self.updatedstats, index=False)


if __name__ == "__main__":
    FilterGff.main()
