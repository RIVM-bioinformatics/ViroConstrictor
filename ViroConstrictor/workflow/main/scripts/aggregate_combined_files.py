import os
from argparse import ArgumentParser
from pathlib import Path

import pandas as pd
from helpers.base_script_class import BaseScript


class AggregateCombinedFiles(BaseScript):
    """
    Aggregate multiple already-processed combined files into a single output file.

    This script concatenates files that have already been processed by earlier
    combination rules (e.g., combine_*_by_sample), which means headers and
    metadata columns are already present.
    """

    def __init__(
        self,
        input: str,
        output: Path | str,
        input_files: list[Path | str],
        file_type: str,
        separator: str = "\t",
    ) -> None:
        super().__init__(input, output)
        self.input_files = input_files
        self.file_type = file_type
        self.separator = separator

    @classmethod
    def add_arguments(cls, parser: ArgumentParser) -> None:
        super().add_arguments(parser)
        parser.add_argument(
            "--input_files",
            nargs="+",
            help="List of input files to aggregate.",
            type=str,
            required=True,
        )
        parser.add_argument(
            "--file_type",
            choices=["mutations", "coverage", "amplicon_coverage"],
            help="Type of file being aggregated.",
            type=str,
            required=True,
        )
        parser.add_argument(
            "--separator",
            help="Field separator (default: tab).",
            type=str,
            default="\t",
        )

    def run(self) -> None:
        """
        Run the aggregation process for combined result files.

        Calls the internal _aggregate_files() method to collect, merge and write
        combined result files into the configured output location.

        Returns
        -------
        None
            This method performs aggregation as a side effect and does not return a value.

        Notes
        -----
        This is a thin public wrapper around the private _aggregate_files implementation.
        """
        self._aggregate_files()

    def _read_files(self, infile_paths: list[Path | str]) -> list[pd.DataFrame]:
        """
        Read multiple CSV files into pandas DataFrames.

        Parameters
        ----------
        infile_paths : list[pathlib.Path or str]
            List of file paths to read. Paths that do not exist or that have size 0 are skipped.

        Returns
        -------
        list[pd.DataFrame]
            List of non-empty DataFrames read from the provided files. Files that are empty
            or raise pandas.errors.EmptyDataError are ignored.

        Raises
        ------
        Exception
            Other exceptions (e.g., file permission errors, parser errors) may propagate.
            pandas.errors.EmptyDataError is caught internally and the corresponding file is skipped.

        Notes
        -----
        - If self.file_type == "amplicon_coverage", files are read with pandas.read_csv(infile).
            Otherwise, files are read with pandas.read_csv(infile, sep=self.separator).
        - The method relies on the instance attributes self.file_type and self.separator.
        """
        dfs = []
        for infile in infile_paths:
            if os.path.exists(infile) and os.path.getsize(infile) > 0:
                try:
                    if self.file_type == "amplicon_coverage":
                        df = pd.read_csv(infile)
                    else:
                        df = pd.read_csv(infile, sep=self.separator)

                    if not df.empty:
                        dfs.append(df)
                except pd.errors.EmptyDataError:

                    continue
        return dfs

    def _aggregate_files(self) -> None:
        """
        Aggregates and concatenates processed input files into a single output file.

        This method reads a list of input files, concatenates their contents (including headers),
        and writes the combined result to the specified output file. If no input files are found,
        it writes an empty file with the appropriate header. The output file format and separator
        depend on the file type.

        Parameters
        ----------
        None

        Returns
        -------
        None

        Notes
        -----
        - For files of type "amplicon_coverage", the output is written as a CSV with default separator.
        - For other file types, the output is written using the specified separator.
        - If no input files are present, an empty file with the correct header is created.
        """
        if dfs := self._read_files(self.input_files):
            combined_df = pd.concat(dfs, ignore_index=True)
            if self.file_type == "amplicon_coverage":
                combined_df.to_csv(self.output, index=False)
            else:
                combined_df.to_csv(self.output, sep=self.separator, index=False)
        else:
            # Write empty file with appropriate header
            self._write_empty_file()

    def _write_empty_file(self) -> None:
        """
        Write an empty output file with appropriate headers based on the file type.

        The method determines the correct header for the file based on the `file_type`
        attribute and writes it to the output file. If the file type is 'amplicon_coverage',
        an empty file is created. For other supported file types, the corresponding header
        line is written.

        Parameters
        ----------
        None

        Returns
        -------
        None

        Notes
        -----
        The headers are predefined for each supported file type:
            - 'mutations': Tab-separated columns for mutation data.
            - 'coverage': Tab-separated columns for coverage statistics.
            - 'amplicon_coverage': No header, creates an empty CSV file.
        """
        headers = {
            "mutations": "Sample\tVirus\tReference_ID\tPosition\tReference_Base\tVariant_Base\tDepth\n",
            "coverage": "Sample_name\tVirus\tReference_ID\tWidth_at_mincov_1\tWidth_at_mincov_5\tWidth_at_mincov_10\tWidth_at_mincov_50\tWidth_at_mincov_100\n",
            "amplicon_coverage": "",  # Empty CSV
        }

        with open(self.output, "w") as f:
            if self.file_type in headers:
                f.write(headers[self.file_type])


if __name__ == "__main__":
    AggregateCombinedFiles.main()
