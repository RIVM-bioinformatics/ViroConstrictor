from argparse import ArgumentParser
from pathlib import Path

from Bio import SeqIO
from helpers.base_script_class import BaseScript  # type: ignore[import]  # noqa: F401,E402


class FilterReferences(BaseScript):
    """
    Filters reference sequences from a FASTA file based on a wildcard segment.

    Parameters
    ----------
    input : str
        Path to the input FASTA file.
    output : str
        Path to the output FASTA file containing the filtered sequences.
    wildcard_segment : str
        Wildcard segment to filter sequences. If "None", all sequences are kept.

    Methods
    -------
    run()
        Executes the filtering of reference sequences.
    """

    def __init__(self, input: Path | str, output: Path | str, wildcard_segment: str) -> None:
        super().__init__(input, output)
        self.wildcard_segment = wildcard_segment

    @classmethod
    def add_arguments(cls, parser: ArgumentParser) -> None:
        super().add_arguments(parser)
        parser.add_argument(
            "--wildcard_segment",
            metavar="String",
            help="Wildcard segment to filter sequences. If 'None', all sequences are kept.",
            type=str,
            required=True,
        )

    def run(self) -> None:
        self._filter_references()

    def _filter_references(self) -> None:
        """
        Filters reference sequences from the input FASTA file based on the wildcard segment
        and writes the filtered sequences to the output FASTA file.
        """
        assert isinstance(self.input, (Path, str)), "Input should be a string path to the FASTA file."
        assert isinstance(self.output, (Path, str)), "Output should be a string path for the filtered FASTA file."
        assert isinstance(self.wildcard_segment, str), "Wildcard segment should be a string."

        # Parse the input FASTA file and filter sequences
        records_to_keep = []
        for record in SeqIO.parse(self.input, "fasta"):
            if self.wildcard_segment != "None":
                if record.description.split(" ")[1].split("|")[0] == self.wildcard_segment:
                    records_to_keep.append(record)
            else:
                records_to_keep.append(record)

        # Write the filtered sequences to the output FASTA file
        SeqIO.write(records_to_keep, self.output, "fasta")


if __name__ == "__main__":
    FilterReferences.main()
