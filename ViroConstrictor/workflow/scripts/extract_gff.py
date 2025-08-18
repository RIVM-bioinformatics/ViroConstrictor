from argparse import ArgumentParser
from pathlib import Path

import AminoExtract

from ViroConstrictor.workflow.scripts import BaseScript


class ExtractGff(BaseScript):
    """
    Extracts a specific GFF record based on the provided reference ID.

    Parameters
    ----------
    input : str
        Path to the input GFF file.
    output : str
        Path to the output file where the extracted GFF record will be saved.
    refID : str
        Reference ID of the GFF record to extract.

    Attributes
    ----------
    input : str
        Input GFF file path.
    output : str
        Output file path for the extracted GFF record.
    refID : str
        Reference ID of the GFF record to extract.
    """

    def __init__(self, input: Path | str, output: Path | str, ref_id: str) -> None:
        super().__init__(input, output)
        self.ref_id = ref_id

    @classmethod
    def add_arguments(cls, parser: ArgumentParser) -> None:
        super().add_arguments(parser)
        parser.add_argument(
            "--refID",
            metavar="String",
            help="Reference ID of the GFF record to extract.",
            type=str,
            required=True,
        )

    def run(self) -> None:
        self._extract_gff_record()

    def _extract_gff_record(self) -> None:
        """
        Extracts a specific GFF record based on the provided reference ID and writes it to the output file.
        """
        assert isinstance(self.input, (Path, str)), "Input should be a string path to the GFF file."
        assert isinstance(self.output, (Path, str)), "Output should be a string path for the extracted GFF record."
        assert isinstance(self.ref_id, str), "Reference ID should be a string."

        gff = AminoExtract.read_gff(self.input)
        gff.df = gff.df[gff.df.seqid == self.ref_id]

        with open(self.output, "w") as f:
            f.write(gff.header)
            f.write(gff.df.to_csv(sep="\t", index=False, header=None))


if __name__ == "__main__":
    ExtractGff.main()
