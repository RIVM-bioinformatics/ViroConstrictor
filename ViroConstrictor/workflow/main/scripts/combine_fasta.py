import os
from argparse import ArgumentParser
from pathlib import Path

from Bio import SeqIO
from helpers.base_script_class import BaseScript


class CombineFasta(BaseScript):
    """Combine FASTA files and modify headers to include Virus and RefID information."""

    def __init__(
        self,
        input: str,
        output: Path | str,
        input_files: list[Path | str],
        virus_list: list[str],
        refid_list: list[str],
    ) -> None:
        super().__init__(input, output)
        self.input_files = input_files
        self.virus_list = virus_list
        self.refid_list = refid_list

    @classmethod
    def add_arguments(cls, parser: ArgumentParser) -> None:
        super().add_arguments(parser)
        parser.add_argument(
            "--input_files",
            nargs="+",
            help="List of input FASTA files to combine.",
            type=str,
            required=True,
        )
        parser.add_argument(
            "--virus_list",
            nargs="+",
            help="List of virus names corresponding to input files.",
            type=str,
            required=True,
        )
        parser.add_argument(
            "--refid_list",
            nargs="+",
            help="List of reference IDs corresponding to input files.",
            type=str,
            required=True,
        )

    def run(self) -> None:
        """
        Executes the process of combining FASTA files.

        This method calls the internal `_combine_fasta` function to perform the
        combination of multiple FASTA files into a single output.

        Returns
        -------
        None
            This method does not return any value.
        """
        self._combine_fasta()

    def _combine_fasta(self) -> None:
        """
        Combine multiple FASTA files into a single output file with modified headers.

        For each input FASTA file, this method reads the sequences and rewrites their headers
        to include the corresponding virus and reference ID information. The new header format is:
        '>sampleID Virus RefID mincov=X' if the original header contains a 'mincov=' field,
        otherwise it appends the virus and reference ID to the original description.

        The combined sequences are written to the specified output file.

        Returns
        -------
        None

        Notes
        -----
        - Skips input files that do not exist or are empty.
        - Assumes that `self.input_files`, `self.virus_list`, and `self.refid_list` are aligned lists.
        - Uses Biopython's SeqIO for FASTA parsing and writing.
        """

        # Create mapping of file to (virus, refid)
        file_map = dict(zip(self.input_files, zip(self.virus_list, self.refid_list)))

        with open(self.output, "w") as out_handle:
            for infile in self.input_files:
                if not os.path.exists(infile) or os.path.getsize(infile) == 0:
                    continue

                virus, refid = file_map[infile]

                for record in SeqIO.parse(infile, "fasta"):
                    # Original header format: >sampleID mincov=X
                    # New format: >sampleID Virus RefID mincov=X
                    header_parts = record.description.split()
                    sample_id = header_parts[0] if header_parts else record.id
                    extra_parts = header_parts[1:] if len(header_parts) > 1 else []
                    if extra_parts and "mincov=" in extra_parts[-1]:
                        mincov_part = extra_parts[-1]
                        # description body should not repeat the sample_id because SeqIO writes ">id description"
                        description_body = f"{virus} {refid} {mincov_part}"
                    elif extra_parts:
                        extra = " ".join(extra_parts)
                        description_body = f"{extra} {virus} {refid}"
                    else:
                        description_body = f"{virus} {refid}"
                    record.id = sample_id
                    record.description = description_body
                    SeqIO.write(record, out_handle, "fasta")


if __name__ == "__main__":
    CombineFasta.main()
