from argparse import ArgumentParser
from pathlib import Path

import pandas as pd

from ViroConstrictor.workflow.scripts.base_script_class import BaseScript


class VcfToTsv(BaseScript):
    def __init__(self, input_path: Path, output_path: Path, samplename: str) -> None:
        super().__init__(input_path, output_path)
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
            self.input_path,
            sep="\t",
            comment="#",
            names=["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"],
        )

        df = df[df.ALT != "N"]

        if df.empty:
            df.to_csv(self.output_path, sep="\t", index=False, header=False)
            return

        df["INFO"] = (
            df["INFO"].str.split("=", expand=True)[1].str.split(";", expand=True)[0]
        )
        df.drop(["ID", "QUAL", "FILTER"], axis=1, inplace=True)
        df.insert(loc=0, column="Sample", value=self.samplename)

        df.to_csv(self.output_path, sep="\t", index=False, header=False)


if __name__ == "__main__":
    VcfToTsv.main()
