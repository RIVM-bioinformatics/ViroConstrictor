"""
Module for grouping amino acid sequences based on sample and feature information.

This module provides functionality to process amino acid sequences from input files,
group them based on metadata from a space file, and write the grouped sequences
to structured output files.
"""

import logging
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

    def __init__(self, input: str, output: str, space: str, log_level: str = "DEBUG") -> None:
        super().__init__(input, output, log_level)
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
        self.log("Starting amino acid grouping process", logging.INFO)

        # Validate input and output arguments
        assert isinstance(self.input, str), "Input should be a space-separated string of file paths."
        assert isinstance(self.output, str), "Output should be a space-separated string of file paths."

        input_files = self.input.split()
        output_files = self.output.split()

        self.log(f"Processing {len(input_files)} input files", logging.INFO)
        self.log(f"Expected {len(output_files)} output files", logging.DEBUG)
        self.log(f"Reading space data from: {self.space}", logging.DEBUG)

        space_data: pd.DataFrame = pd.read_pickle(self.space)
        self.log(f"Loaded space data with {len(space_data)} records", logging.INFO)

        def feature_in_id(feature: str, record_id: str) -> bool:
            pattern = rf"(?:^|[.\-_ ]){re.escape(feature)}(?:$|[.\-_ ])"
            return re.search(pattern, record_id, re.IGNORECASE) is not None

        def process_record(record: SeqIO.SeqRecord, space_data: pd.DataFrame) -> pd.DataFrame:
            sample = record.id.split(".")[0]
            self.log(f"Processing record: {record.id} (sample: {sample})", logging.DEBUG)

            sample_data = space_data.loc[space_data["sample"] == sample]
            if sample_data.empty:
                self.log(f"No sample data found for sample: {sample}", logging.WARNING)
                return pd.DataFrame()

            sample_data = sample_data.explode("AA_FEAT_NAMES").reset_index(drop=True)
            sample_data = sample_data[["sample", "Virus", "RefID", "AA_FEAT_NAMES"]]

            # this filter is wrong because there are many ORF's which are features that are not in the seq_record id.
            all_unique_features = list(sample_data["AA_FEAT_NAMES"].unique())
            self.log(f"Found {len(all_unique_features)} unique features for sample {sample}", logging.DEBUG)

            for feature in all_unique_features:
                if feature_in_id(feature, record.id):
                    self.log(f"Matched feature '{feature}' in record ID: {record.id}", logging.DEBUG)
                    return sample_data.loc[sample_data["AA_FEAT_NAMES"] == feature].assign(AA_SEQ=str(record.seq))

            aa_feature = ".".join(record.id.split(".")[1:])
            self.log(f"Using fallback feature matching: {aa_feature}", logging.DEBUG)
            result = sample_data.loc[sample_data["AA_FEAT_NAMES"] == aa_feature].assign(AA_SEQ=str(record.seq))

            if result.empty:
                self.log(f"No feature match found for record: {record.id}", logging.WARNING)

            return result

        def write_sequences(
            virus: str,
            ref_id: str,
            feature_name: str,
            output_files: list[str],
            seq_records_df: pd.DataFrame,
        ) -> None:
            output_path = f"results/Virus~{virus}/RefID~{ref_id}/aminoacids/{feature_name}.faa"

            self.log(f"Attempting to write sequences for {virus}/{ref_id}/{feature_name}", logging.DEBUG)
            self.log(f"Constructed output path: {output_path}", logging.DEBUG)
            self.log(f"Output files list has {len(output_files)} entries", logging.DEBUG)

            # Debug: Show first few expected output files for comparison
            if output_files:
                self.log(f"First few expected output files: {output_files[:3]}", logging.DEBUG)

            if output_path in output_files:
                filtered_data = seq_records_df.loc[
                    (seq_records_df["Virus"] == virus) & (seq_records_df["RefID"] == ref_id) & (seq_records_df["AA_FEAT_NAMES"] == feature_name)
                ]

                self.log(f"Writing {len(filtered_data)} sequences to: {output_path}", logging.INFO)

                if filtered_data.empty:
                    self.log(f"No sequences found for {virus}/{ref_id}/{feature_name}", logging.WARNING)
                    return

                file_path = Path(output_path)
                self.log(f"Creating directory: {file_path.parent}", logging.DEBUG)
                file_path.parent.mkdir(parents=True, exist_ok=True)

                try:
                    with open(file_path, "w", encoding="utf-8") as file:
                        for _, row in filtered_data.iterrows():
                            file.write(f">{row['sample']}.{row['AA_FEAT_NAMES']}\n{row['AA_SEQ']}\n")
                    self.log(f"Successfully wrote sequences to: {output_path}", logging.DEBUG)

                    # Verify file was actually created
                    if file_path.exists():
                        file_size = file_path.stat().st_size
                        self.log(f"File verification: {output_path} exists with size {file_size} bytes", logging.DEBUG)
                    else:
                        self.log(f"File verification FAILED: {output_path} does not exist after write!", logging.ERROR)

                except Exception as e:
                    self.log(f"Error writing to file {output_path}: {str(e)}", logging.ERROR)
                    raise
            else:
                self.log(f"Output path not in expected files: {output_path}", logging.WARNING)
                self.log(f"Looking for exact matches in output_files list...", logging.DEBUG)

                # Debug: Look for similar paths
                similar_paths = [f for f in output_files if feature_name in f and virus in f and ref_id in f]
                if similar_paths:
                    self.log(f"Found similar paths: {similar_paths}", logging.DEBUG)
                else:
                    self.log(f"No similar paths found containing {feature_name}, {virus}, {ref_id}", logging.DEBUG)

        seq_records_df = pd.DataFrame()
        total_records_processed = 0

        for input_file in input_files:
            self.log(f"Processing input file: {input_file}", logging.INFO)
            try:
                records: list[SeqIO.SeqRecord] = list(SeqIO.parse(input_file, "fasta"))
                self.log(f"Found {len(records)} records in {input_file}", logging.DEBUG)

                for record in records:
                    processed_data = process_record(record, space_data)
                    if not processed_data.empty:
                        seq_records_df = pd.concat([seq_records_df, processed_data], ignore_index=True)
                        total_records_processed += 1

            except Exception as e:
                self.log(f"Error processing file {input_file}: {str(e)}", logging.ERROR)
                raise

        self.log(f"Processed {total_records_processed} records total", logging.INFO)

        if seq_records_df.empty:
            self.log("No valid sequence records found after processing", logging.WARNING)
            return

        unique_viruses = np.array(seq_records_df["Virus"].unique(), dtype=np.str_)
        self.log(f"Found {len(unique_viruses)} unique viruses: {', '.join(unique_viruses)}", logging.INFO)

        files_written = 0
        for virus in unique_viruses:
            unique_ref_ids = np.array(
                seq_records_df.loc[seq_records_df["Virus"] == virus]["RefID"].unique(),
                dtype=np.str_,
            )
            self.log(f"Processing virus '{virus}' with {len(unique_ref_ids)} reference IDs", logging.DEBUG)

            for ref_id in unique_ref_ids:
                unique_feature_names = np.array(
                    seq_records_df.loc[(seq_records_df["Virus"] == virus) & (seq_records_df["RefID"] == ref_id)]["AA_FEAT_NAMES"].unique(),
                    dtype=np.str_,
                )

                for feature_name in unique_feature_names:
                    write_sequences(virus, ref_id, feature_name, output_files, seq_records_df)
                    files_written += 1

        self.log(f"Amino acid grouping completed successfully. Wrote {files_written} output files.", logging.INFO)


if __name__ == "__main__":
    GroupAminoAcids.main()
