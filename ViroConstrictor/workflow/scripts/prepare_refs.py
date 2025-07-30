from argparse import ArgumentParser
from pathlib import Path

from Bio import SeqIO

from ViroConstrictor.workflow.scripts.base_script_class import BaseScript


class PrepareRefs(BaseScript):
    """
    PrepareRefs script to process reference sequences.

    This script reads a FASTA file, converts the sequence to uppercase,
    and writes it to an output file if the reference ID matches.
    """

    def __init__(self, input_path: Path, output_path: Path, reference_id: str) -> None:
        super().__init__(input_path, output_path)
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
        for record in SeqIO.parse(self.input_path.as_posix(), "fasta"):
            if self.reference_id in record.id:
                record.seq = record.seq.upper()
                SeqIO.write(record, self.output_path.as_posix(), "fasta")


if __name__ == "__main__":
    PrepareRefs.main()
