from argparse import ArgumentParser
from pathlib import Path
from typing import Iterator, cast

from Bio import SeqIO

from .base_script_class import BaseScript


class PrepareRefs(BaseScript):
    """
    PrepareRefs script to process reference sequences.

    This script reads a FASTA file, converts the sequence to uppercase,
    and writes it to an output file if the reference ID matches.
    """

    def __init__(self, input: Path, output: Path, reference_id: str) -> None:
        print("HOIASDS")
        super().__init__(input, output)

        self.reference_id = reference_id

    @classmethod
    def add_arguments(cls, parser: ArgumentParser) -> None:
        super().add_arguments(parser)
        parser.add_argument(
            "--reference_id",
            metavar="String",
            help="Reference ID to match in the FASTA file.",
            type=str,
            required=True,
        )

    def run(self) -> None:
        self._prepare_refs()

    def _prepare_refs(self) -> None:
        """Prepare reference sequences by converting to uppercase."""
        assert not isinstance(self.input, list), "Input must be cannot be a list of strs."
        assert not isinstance(self.output, list), "Output must be cannot be a list of strs."
        records = cast(Iterator[SeqIO.SeqRecord], SeqIO.parse(self.input, "fasta"))
        for record in records:
            if self.reference_id in record.id:
                record.seq = record.seq.upper()
                SeqIO.write(record, self.output, "fasta")


if __name__ == "__main__":
    PrepareRefs.main()
