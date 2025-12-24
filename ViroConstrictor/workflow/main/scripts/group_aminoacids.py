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
        """
        Adds command-line arguments specific to this script.
        First the common arguments from BaseScript are added.
        Then the specific argument for this script is added.
        """
        super().add_arguments(parser)
        parser.add_argument(
            "--space",
            metavar="List of Files",
            help="Path to the space file containing sample and amino acid feature names.",
            type=Path,
            required=True,
        )

    def run(self) -> None:
        """Executes the grouping of amino acid sequences."""
        self._group_amino_acids()

    def _debug_space_data(self, space_data: pd.DataFrame) -> None:
        """Logs debug information about the space data structure."""
        if not space_data.empty:
            self.log(logging.DEBUG, "Space data columns: %s", list(space_data.columns))
            # Show unique AA_FEAT_NAMES to see what features are expected
            if "AA_FEAT_NAMES" in space_data.columns:
                all_features: list[str] = []
                for _, row in space_data.iterrows():
                    if pd.notna(row["AA_FEAT_NAMES"]):
                        if isinstance(row["AA_FEAT_NAMES"], list):
                            all_features.extend(row["AA_FEAT_NAMES"])
                        else:
                            all_features.append(row["AA_FEAT_NAMES"])
                unique_features = list(set(all_features))
                self.log(logging.DEBUG, "Unique AA_FEAT_NAMES from space data: %s", unique_features)

            # Show sample of space data
            self.log(logging.DEBUG, "Sample space data:\n%s", space_data.head())

    def _group_amino_acids(self) -> None:
        """
        Groups amino acid sequences from input files based on sample and feature information
        and writes them to output files in a structured directory format.
        """
        self.log(logging.INFO, "Starting amino acid grouping process")

        # Validate input and output arguments
        assert isinstance(self.input, str), "Input should be a space-separated string of file paths."
        assert isinstance(self.output, str), "Output should be a space-separated string of file paths."

        input_files = self.input.split()
        output_files = self.output.split()

        self.log(logging.INFO, "Processing %s input files", len(input_files))
        self.log(logging.DEBUG, "Expected %s output files", len(output_files))
        self.log(logging.DEBUG, "Reading space data from: %s", self.space)

        space_data: pd.DataFrame = pd.read_pickle(self.space)
        self.log(logging.INFO, "Loaded space data with %s records", len(space_data))

        if not space_data.empty:
            self.log(logging.DEBUG, "Space data columns: %s", list(space_data.columns))

            # Show sample of space data
            self.log(logging.DEBUG, "Sample space data:\n%s", space_data.head())

        def feature_in_id(feature: str, record_id: str) -> bool:
            pattern = rf"(?:^|[.\-_ ]){re.escape(feature)}(?:$|[.\-_ ])"
            return re.search(pattern, record_id, re.IGNORECASE) is not None

        def process_record(record: SeqIO.SeqRecord, space_data: pd.DataFrame) -> pd.DataFrame:
            sample = record.id.split(".")[0]
            self.log(logging.DEBUG, "Processing record: %s (sample: %s)", record.id, sample)

            sample_data = space_data.loc[space_data["sample"] == sample]
            if sample_data.empty:
                self.log(logging.WARNING, "No sample data found for sample: %s", sample)
                return pd.DataFrame()

            sample_data = sample_data.explode("AA_FEAT_NAMES").reset_index(drop=True)
            sample_data = sample_data[["sample", "Virus", "RefID", "AA_FEAT_NAMES"]]

            # this filter is wrong because there are many ORF's which are features that are not in the seq_record id.
            all_unique_features = list(sample_data["AA_FEAT_NAMES"].unique())
            self.log(logging.DEBUG, f"Found {len(all_unique_features)} unique features for sample {sample}")

            for feature in all_unique_features:
                if feature_in_id(feature, record.id):
                    self.log(logging.DEBUG, "Matched feature '%s' in record ID: %s", feature, record.id)
                    return sample_data.loc[sample_data["AA_FEAT_NAMES"] == feature].assign(AA_SEQ=str(record.seq))

            aa_feature = ".".join(record.id.split(".")[1:])
            self.log(logging.DEBUG, "Using fallback feature matching: %s", aa_feature)
            result = sample_data.loc[sample_data["AA_FEAT_NAMES"] == aa_feature].assign(AA_SEQ=str(record.seq))

            if result.empty:
                self.log(logging.WARNING, "No feature match found for record: %s", record.id)
            return result

        def write_sequences(
            virus: str,
            ref_id: str,
            feature_name: str,
            output_files: list[str],
            seq_records_df: pd.DataFrame,
        ) -> None:
            output_path = f"results/Virus~{virus}/RefID~{ref_id}/aminoacids/{feature_name}.faa"

            self.log(logging.DEBUG, "Attempting to write sequences for %s/%s/%s", virus, ref_id, feature_name)
            self.log(logging.DEBUG, "Constructed output path: %s", output_path)
            self.log(logging.DEBUG, "Output files list has %d entries", len(output_files))

            # Debug: Show first few expected output files for comparison
            if output_files:
                self.log(logging.DEBUG, "First few expected output files: %s", output_files[:3])
            if output_path in output_files:
                filtered_data = seq_records_df.loc[
                    (seq_records_df["Virus"] == virus) & (seq_records_df["RefID"] == ref_id) & (seq_records_df["AA_FEAT_NAMES"] == feature_name)
                ]

                self.log(logging.INFO, "Writing %d sequences to: %s", len(filtered_data), output_path)

                if filtered_data.empty:
                    self.log(logging.WARNING, "No sequences found for %s/%s/%s", virus, ref_id, feature_name)
                    return

                file_path = Path(output_path)
                self.log(logging.DEBUG, "Creating directory: %s", file_path.parent)
                file_path.parent.mkdir(parents=True, exist_ok=True)

                try:
                    with open(file_path, "w", encoding="utf-8") as file:
                        for _, row in filtered_data.iterrows():
                            file.write(f">{row['sample']}.{row['AA_FEAT_NAMES']}\n{row['AA_SEQ']}\n")
                    self.log(logging.DEBUG, "Successfully wrote sequences to: %s", output_path)

                    # Verify file was actually created
                    if file_path.exists():
                        file_size = file_path.stat().st_size
                        self.log(logging.DEBUG, "File verification: %s exists with size %d bytes", output_path, file_size)
                    else:
                        self.log(logging.ERROR, "File verification FAILED: %s does not exist after write!", output_path)

                except Exception as e:
                    self.log(logging.ERROR, "Error writing to file %s: %s", output_path, str(e))
                    raise
            else:
                self.log(logging.WARNING, "Output path not in expected files: %s", output_path)
                self.log(logging.DEBUG, "Looking for exact matches in output_files list...")
                # Debug: Look for similar paths
                similar_paths = [f for f in output_files if feature_name in f and virus in f and ref_id in f]
                if similar_paths:
                    self.log(logging.DEBUG, "Found similar paths: %s", similar_paths)
                else:
                    self.log(logging.DEBUG, "No similar paths found containing %s, %s, %s", feature_name, virus, ref_id)

        seq_records_df = pd.DataFrame()
        total_records_processed = 0

        for input_file in input_files:
            self.log(logging.INFO, "Processing input file: %s", input_file)
            try:
                records: list[SeqIO.SeqRecord] = list(SeqIO.parse(input_file, "fasta"))
                self.log(logging.DEBUG, "Found %d records in %s", len(records), input_file)

                for record in records:
                    processed_data = process_record(record, space_data)
                    if not processed_data.empty:
                        seq_records_df = pd.concat([seq_records_df, processed_data], ignore_index=True)
                        total_records_processed += 1

            except Exception as e:
                self.log(logging.ERROR, "Error processing file %s: %s", input_file, str(e))
                raise

        self.log(logging.INFO, "Processed %d records total", total_records_processed)

        if seq_records_df.empty:
            self.log(logging.WARNING, "No valid sequence records found after processing")
            return

        unique_viruses = np.array(seq_records_df["Virus"].unique(), dtype=np.str_)
        self.log(logging.INFO, "Found %d unique viruses: %s", len(unique_viruses), ", ".join(unique_viruses))
        files_written = 0
        total_expected_files = len(output_files)
        files_written_details = []

        for virus in unique_viruses:
            unique_ref_ids = np.array(
                seq_records_df.loc[seq_records_df["Virus"] == virus]["RefID"].unique(),
                dtype=np.str_,
            )
            self.log(logging.DEBUG, "Processing virus '%s' with %d reference IDs", virus, len(unique_ref_ids))
            self.log(logging.DEBUG, "Reference IDs for virus %s: %s", virus, list(unique_ref_ids))

            for ref_id in unique_ref_ids:
                unique_feature_names = np.array(
                    seq_records_df.loc[(seq_records_df["Virus"] == virus) & (seq_records_df["RefID"] == ref_id)]["AA_FEAT_NAMES"].unique(),
                    dtype=np.str_,
                )
                self.log(logging.DEBUG, "Processing RefID '%s' with features: %s", ref_id, list(unique_feature_names))

                for feature_name in unique_feature_names:
                    expected_path = f"results/Virus~{virus}/RefID~{ref_id}/aminoacids/{feature_name}.faa"
                    self.log(logging.DEBUG, "About to write feature '%s' to expected path: %s", feature_name, expected_path)

                    # Check if this path is in the expected output files
                    path_in_expected = expected_path in output_files
                    self.log(logging.DEBUG, "Expected path '%s' is in output_files: %s", expected_path, path_in_expected)

                    write_sequences(virus, ref_id, feature_name, output_files, seq_records_df)
                    files_written += 1
                    files_written_details.append(f"{virus}/{ref_id}/{feature_name}")

        self.log(logging.DEBUG, "Files written details: %s", files_written_details)
        self.log(logging.INFO, "Total files written: %d, Total expected: %d", files_written, total_expected_files)
        self.log(logging.INFO, "Amino acid grouping completed successfully. Wrote %d output files.", files_written)


if __name__ == "__main__":
    GroupAminoAcids.main()
