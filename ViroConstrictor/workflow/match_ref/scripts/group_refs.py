from argparse import ArgumentParser
from pathlib import Path

import pandas as pd
from Bio import SeqIO

from ViroConstrictor.workflow.helpers.base_script_class import BaseScript  # type: ignore[import]  # noqa: F401,E402


class GroupRefs(BaseScript):
    """
    Groups reference sequences and their associated statistics into a single FASTA file and CSV file.

    Parameters
    ----------
    input_refs : list[str]
        List of input FASTA file paths containing reference sequences.
    input_stats : list[str]
        List of input CSV file paths containing reference statistics.
    output_ref : str
        Path to the output FASTA file containing grouped reference sequences.
    output_stats : str
        Path to the output CSV file containing grouped reference statistics.
    sample : str
        Sample name to associate with the grouped references.

    Methods
    -------
    run()
        Executes the grouping of reference sequences and statistics.
    """

    def __init__(
        self,
        input: str,
        input_refs: list[Path | str],  # this might be a wrong type assignment
        input_stats: list[Path | str],
        output: Path | str,
        output_stats: Path | str,
        sample: str,
    ) -> None:
        super().__init__(input_refs, output)
        self.input_stats = input_stats
        self.output_stats = output_stats
        self.sample = sample

    @classmethod
    def add_arguments(cls, parser: ArgumentParser) -> None:
        super().add_arguments(parser)
        parser.add_argument(
            "--input_refs",
            metavar="File",
            nargs="+",
            help="List of input FASTA file paths containing reference sequences.",
            required=True,
        )
        parser.add_argument(
            "--input_stats",
            metavar="File",
            nargs="+",
            help="List of input CSV file paths containing reference statistics.",
            required=True,
        )
        parser.add_argument(
            "--output_stats",
            metavar="File",
            help="Path to the output CSV file containing grouped reference statistics.",
            required=True,
        )
        parser.add_argument(
            "--sample",
            metavar="String",
            help="Sample name to associate with the grouped references.",
            required=True,
        )

    def run(self) -> None:
        self._group_refs()

    def _group_refs(self) -> None:
        """
        Groups reference sequences and their associated statistics into a single FASTA file and CSV file.
        """
        assert isinstance(self.input, list), "Input_refs should be a list of file paths."
        assert isinstance(self.input_stats, list), "Input_stats should be a list of file paths."
        assert isinstance(self.output, (Path, str)), "Output_ref should be a string path for the output FASTA file."
        assert isinstance(self.output_stats, (Path, str)), "Output_stats should be a string path for the output CSV file."
        assert isinstance(self.sample, str), "Sample should be a string."

        # Read and combine reference sequences
        seqrecords = []
        for file in self.input:
            seqrecords.extend(iter(SeqIO.parse(file, "fasta")))

        # Read and combine statistics
        df = pd.DataFrame()
        for file in self.input_stats:
            tempdf = pd.read_csv(file)
            df = pd.concat([df, tempdf], ignore_index=True)

        df["sample"] = self.sample

        # Rename sequence records for segmented or non-segmented mode
        renamed_seqrecords = []
        if len(seqrecords) > 1:
            # Segmented mode
            for record in seqrecords:
                record.id = record.description.split()[1].split("|")[0]
                record.description = " ".join([record.name, " ".join(record.description.split(" ")[1:])])
                renamed_seqrecords.append(record)
        else:
            # Non-segmented mode
            renamed_seqrecords.append(seqrecords[0])

        # Add sequence record information to the dataframe
        for record in renamed_seqrecords:
            df.loc[df["Reference"].str.contains(record.id), "seqrecord_id"] = record.id
            df.loc[df["Reference"].str.contains(record.id), "seqrecord_description"] = record.description
            df.loc[df["Reference"].str.contains(record.id), "seqrecord_name"] = record.name
            df.loc[df["Reference"].str.contains(record.id), "seqrecord_seq"] = str(record.seq)

        # Replace the "Reference" column with the "seqrecord_id" column
        df["Reference"] = df["seqrecord_id"]
        df["Reference_file"] = Path(self.output).resolve()
        # Write the dataframe to a CSV file
        df.to_csv(self.output_stats, index=False)

        # Write the renamed sequence records to the output FASTA file
        SeqIO.write(renamed_seqrecords, self.output, "fasta")


if __name__ == "__main__":
    GroupRefs.main()
