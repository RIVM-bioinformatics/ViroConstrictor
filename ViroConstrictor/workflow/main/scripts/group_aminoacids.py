"""
Module for grouping amino acid sequences based on sample and feature information.

This module provides functionality to process amino acid sequences from input files,
group them based on metadata from a space file, and write the grouped sequences
to structured output files.
"""

import re
from argparse import ArgumentParser
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO

from helpers.base_script_class import BaseScript  # type: ignore[import]  # noqa: F401,E402


class GroupAminoAcids(BaseScript):
    """
    A script to group amino acid sequences based on sample and feature information.

    This class reads amino acid sequences from input FASTA files, matches them with
    metadata from a space file, and writes the grouped sequences to output files
    in a structured directory format.

    Parameters
    ----------
    input : str
        A space-separated string of input file paths containing amino acid sequences in FASTA format.
    output : str
        A space-separated string of output file paths where grouped sequences will be written.
    space : Path
        Path to the space file containing metadata about samples and amino acid feature names.

    Methods
    -------
    add_arguments(parser: ArgumentParser)
        Adds command-line arguments specific to this script.
    run()
        Executes the grouping of amino acid sequences.
    _group_amino_acids()
        Internal method to perform the grouping and writing of sequences.
    """

    def __init__(self, input: str, output: str, space: str) -> None:
        super().__init__(input, output)
        self.space = Path(space)

    @classmethod
    def add_arguments(cls, parser: ArgumentParser):
        super().add_arguments(parser)
        parser.add_argument(
            "--space",
            metavar="List of Files",
            help="Path to the space file containing sample and amino acid feature names.",
            type=Path,
            required=True,
        )

    def run(self) -> None:
        self._group_amino_acids()

    def _group_amino_acids(self) -> None:
        """
        Groups amino acid sequences from input files based on sample and feature information
        and writes them to output files in a structured directory format.
        """
        # Validate input and output arguments
        assert isinstance(
            self.input, str
        ), "Input should be a space-separated string of file paths."
        assert isinstance(
            self.output, str
        ), "Output should be a space-separated string of file paths."

        input_files = self.input.split()
        output_files = self.output.split()
        space_data: pd.DataFrame = pd.read_pickle(self.space)

        def feature_in_id(feature: str, record_id: str) -> bool:
            pattern = rf"(?:^|[.\-_ ]){re.escape(feature)}(?:$|[.\-_ ])"
            return re.search(pattern, record_id, re.IGNORECASE) is not None

        def process_record(
            record: SeqIO.SeqRecord, space_data: pd.DataFrame
        ) -> pd.DataFrame:
            sample = record.id.split(".")[0]
            sample_data = space_data.loc[space_data["sample"] == sample]
            sample_data = sample_data.explode("AA_FEAT_NAMES").reset_index(drop=True)
            sample_data = sample_data[["sample", "Virus", "RefID", "AA_FEAT_NAMES"]]

            # this filter is wrong because there are many ORF's which are features that are not in the seq_record id.
            all_unique_features = list(sample_data["AA_FEAT_NAMES"].unique())

            for feature in all_unique_features:
                if feature_in_id(feature, record.id):
                    return sample_data.loc[
                        sample_data["AA_FEAT_NAMES"] == feature
                    ].assign(AA_SEQ=str(record.seq))

            aa_feature = ".".join(record.id.split(".")[1:])
            return sample_data.loc[sample_data["AA_FEAT_NAMES"] == aa_feature].assign(
                AA_SEQ=str(record.seq)
            )

        def write_sequences(
            virus: str,
            ref_id: str,
            feature_name: str,
            output_files: list[str],
            seq_records_df: pd.DataFrame,
        ) -> None:
            output_path = (
                f"results/Virus~{virus}/RefID~{ref_id}/aminoacids/{feature_name}.faa"
            )
            if output_path in output_files:
                filtered_data = seq_records_df.loc[
                    (seq_records_df["Virus"] == virus)
                    & (seq_records_df["RefID"] == ref_id)
                    & (seq_records_df["AA_FEAT_NAMES"] == feature_name)
                ]
                file_path = Path(output_path)
                file_path.parent.mkdir(parents=True, exist_ok=True)
                with open(file_path, "w", encoding="utf-8") as file:
                    for _, row in filtered_data.iterrows():
                        file.write(
                            f">{row['sample']}.{row['AA_FEAT_NAMES']}\n{row['AA_SEQ']}\n"
                        )

        seq_records_df = pd.DataFrame()
        for input_file in input_files:
            records: list[SeqIO.SeqRecord] = list(SeqIO.parse(input_file, "fasta"))
            for record in records:
                processed_data = process_record(record, space_data)
                seq_records_df = pd.concat(
                    [seq_records_df, processed_data], ignore_index=True
                )

        unique_viruses = np.array(seq_records_df["Virus"].unique(), dtype=np.str_)
        for virus in unique_viruses:
            unique_ref_ids = np.array(
                seq_records_df.loc[seq_records_df["Virus"] == virus]["RefID"].unique(),
                dtype=np.str_,
            )
            for ref_id in unique_ref_ids:
                unique_feature_names = np.array(
                    seq_records_df.loc[
                        (seq_records_df["Virus"] == virus)
                        & (seq_records_df["RefID"] == ref_id)
                    ]["AA_FEAT_NAMES"].unique(),
                    dtype=np.str_,
                )
                for feature_name in unique_feature_names:

                    write_sequences(
                        virus, ref_id, feature_name, output_files, seq_records_df
                    )


if __name__ == "__main__":
    GroupAminoAcids.main()
