from argparse import ArgumentParser

import pandas as pd
from helpers.base_script_class import BaseScript  # type: ignore[import]  # noqa: F401,E402


class FilterBedInput(BaseScript):
    """
    FilterBedInput script to filter BED file based on reference ID.

    This script reads a BED file, filters it by a specified reference ID,
    and writes the filtered data to an output file.
    """

    def __init__(self, input: str, output: str, reference_id: str, log_level: str = "INFO") -> None:
        super().__init__(input, output, log_level)
        self.reference_id = reference_id

    @classmethod
    def add_arguments(cls, parser: ArgumentParser) -> None:
        super().add_arguments(parser)
        parser.add_argument(
            "--reference_id",
            metavar="String",
            help="Reference ID to filter the BED file.",
            type=str,
            required=True,
        )

    def run(self) -> None:
        self._filter_bed_input()

    def _filter_bed_input(self) -> None:
        """Filter BED input based on reference ID."""
        assert isinstance(self.input, str), "Input must be a string path."
        assert isinstance(self.output, str), "Output must be a string path."

        df = pd.read_csv(
            self.input,
            sep="\t",
            names=["ref", "start", "stop", "name", "score", "strand"],
            index_col=False,
            keep_default_na=False,
        )

        df = df[df.ref == self.reference_id]

        df.to_csv(self.output, sep="\t", index=False, header=False)


if __name__ == "__main__":
    FilterBedInput.main()
