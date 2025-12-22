from argparse import ArgumentParser
from pathlib import Path

import pandas as pd
from helpers.base_script_class import BaseScript  # type: ignore[import]  # noqa: F401,E402


class VcfToTsv(BaseScript):
    def __init__(self, input: Path, output: Path, samplename: str, log_level: str = "INFO") -> None:
        super().__init__(input, output, log_level)
        self.samplename = samplename

    @classmethod
    def add_arguments(cls, parser: ArgumentParser) -> None:
        super().add_arguments(parser)
        parser.add_argument(
            "--samplename",
            metavar="String",
            help="Sample name to be added to the output.",
            type=str,
            required=True,
        )

    def run(self) -> None:
        self._convert_vcf_to_tsv()

    def _convert_vcf_to_tsv(self) -> None:
        """Convert VCF to TSV format."""
        df = pd.read_csv(
            self.input,
            sep="\t",
            comment="#",
            names=["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"],
        )

        df = df[df.ALT != "N"]

        if df.empty:
            df.to_csv(self.output, sep="\t", index=False, header=False)
            return

        df["INFO"] = df["INFO"].str.split("=", expand=True)[1].str.split(";", expand=True)[0]
        df.drop(["ID", "QUAL", "FILTER"], axis=1, inplace=True)
        df.insert(loc=0, column="Sample", value=self.samplename)

        df.to_csv(self.output, sep="\t", index=False, header=False)


if __name__ == "__main__":
    VcfToTsv.main()
