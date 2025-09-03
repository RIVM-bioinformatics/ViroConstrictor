from argparse import ArgumentParser
from pathlib import Path

import pandas as pd
from Bio import SeqIO

from .base_script_class import BaseScript


class Boc(BaseScript):
    def __init__(
        self,
        input: Path | str,
        output: Path | str,
        samplename: str,
        coverage: Path | str,
    ) -> None:
        super().__init__(input, output)
        self.samplename = samplename
        self.coverages = coverage

    @classmethod
    def add_arguments(cls, parser: ArgumentParser) -> None:
        super().add_arguments(parser)
        parser.add_argument(
            "--samplename",
            metavar="String",
            help="Name of the sample.",
            type=str,
            required=True,
        )
        parser.add_argument(
            "--coverage",
            metavar="String",
            help="Path to the coverage file.",
            type=str,
            required=True,
        )

    def run(self) -> None:
        self._boc()

    def _boc(self) -> None:
        assert isinstance(
            self.input, (Path, str)
        ), "Input should be a string path to the GFF file."
        assert isinstance(
            self.output, (Path, str)
        ), "Output should be a string path for the extracted GFF record."

        df = pd.read_csv(self.coverages, sep="\t", index_col=0)

        for record in SeqIO.parse(self.input, "fasta"):
            a = len(record.seq)

        with open(self.output, "w") as f:
            f.write(
                f"{self.samplename}\t{100-(int(df[df < 1].count().iloc[0])/a)*100}\t{100-(int(df[df < 5].count().iloc[0])/a)*100}\t{100-(int(df[df < 10].count().iloc[0])/a)*100}\t{100-(int(df[df < 50].count().iloc[0])/a)*100}\t{100-(int(df[df < 100].count().iloc[0])/a)*100}\n"
            )


if __name__ == "__main__":
    Boc.main()
