import os
from argparse import ArgumentParser
from pathlib import Path
from typing import Generator, Literal

import pandas as pd
from helpers.base_script_class import BaseScript


class CombineTabular(BaseScript):
    """Combine TSV/CSV files and add Virus and RefID metadata columns."""

    def __init__(
        self,
        input: str,
        output: Path | str,
        input_files: list[Path | str],
        virus_list: list[str],
        refid_list: list[str],
        file_type: Literal["mutations", "coverage", "amplicon_coverage"],
        separator: str = "\t",
    ) -> None:
        super().__init__(input, output)
        self.input_files = input_files
        self.virus_list = virus_list
        self.refid_list = refid_list
        self.file_type = file_type
        self.separator = separator

    @classmethod
    def add_arguments(cls, parser: ArgumentParser) -> None:
        super().add_arguments(parser)
        parser.add_argument(
            "--input_files",
            nargs="+",
            help="List of input files to combine.",
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
        parser.add_argument(
            "--file_type",
            choices=["mutations", "coverage", "amplicon_coverage"],
            help="Type of file being combined (determines header format).",
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
        Executes the tabular combination process.

        This method calls the internal `_combine_tabular` function to combine tabular data as part of the workflow.

        Returns
        -------
        None
            This method does not return any value.
        """
        self._combine_tabular()

    def _read_tabular_file(
        self,
        infile: Path | str,
        sep: str | None = None,
        header: int | None = 0,
        names: list[str] | None = None,
    ) -> pd.DataFrame | None:
        """
        Read a tabular file with configurable parameters.

        Parameters
        ----------
        infile : Path | str
            Path to the input file.
        sep : str | None
            Field separator. If None, uses self.separator.
        header : int | None
            Row number to use as column names. Use None for no header.
        names : list[str] | None
            Column names to use when header is None.

        Returns
        -------
        pd.DataFrame | None
            DataFrame if file is non-empty, None otherwise.
        """
        if sep is None:
            sep = self.separator

        try:
            df = pd.read_csv(infile, sep=sep, header=header, names=names)
        except pd.errors.EmptyDataError:
            return None
        return None if df.empty else df

    def _write_empty_tabular(self, fileheader: str) -> None:
        """
        Write an empty tabular file with the appropriate header based on the file type.
        Depending on the `fileheader` argument, this method writes a tab-separated header line
        to the output file specified by `self.output`. If `fileheader` is "coverage", it writes
        the coverage columns. If `fileheader` is "mutations", it writes the mutation columns.
        Parameters
        ----------
        fileheader : str
            The type of tabular file to write. Accepted values are "coverage" and "mutations".
        Returns
        -------
        None
            This method does not return any value. It writes directly to the output file.
        """

        columns = ""
        if fileheader == "coverage":
            columns = (
                "Sample_name\tVirus\tReference_ID\tWidth_at_mincov_1\tWidth_at_mincov_5"
                "\tWidth_at_mincov_10\tWidth_at_mincov_50\tWidth_at_mincov_100\n"
            )
        elif fileheader == "mutations":
            columns = "Sample\tVirus\tReference_ID\tPosition\tReference_Base\tVariant_Base\tDepth\n"

        with open(self.output, "w") as f:
            f.write(columns)

    def _reorder_mutations_columns(self, combined_df: pd.DataFrame) -> pd.DataFrame:
        """
        Reorders the columns of the combined DataFrame so that the 'Virus' column directly follows the 'Sample' column.

        Parameters
        ----------
        combined_df : pd.DataFrame
            The DataFrame containing mutation data with at least 'Sample' and 'Virus' columns.

        Returns
        -------
        pd.DataFrame
            The DataFrame with the 'Virus' column positioned immediately after the 'Sample' column.
            If either 'Sample' or 'Virus' is not present, returns the DataFrame unchanged.
        """
        cols = combined_df.columns.tolist()
        if "Sample" in cols and "Virus" in cols:
            cols.remove("Virus")
            sample_idx = cols.index("Sample")
            cols.insert(sample_idx + 1, "Virus")
            return combined_df[cols]
        return combined_df

    def _reorder_amplicon_columns(self, combined_df: pd.DataFrame) -> pd.DataFrame:
        """
        Reorders the columns of the given DataFrame to move 'Virus' and 'Reference_ID' to the front.

        Parameters
        ----------
        combined_df : pd.DataFrame
            The DataFrame whose columns are to be reordered.

        Returns
        -------
        pd.DataFrame
            A DataFrame with 'Virus' and 'Reference_ID' columns moved to the front if both are present;
            otherwise, returns the original DataFrame unchanged.
        """
        cols = combined_df.columns.tolist()
        if "Virus" in cols and "Reference_ID" in cols:
            cols.remove("Virus")
            cols.remove("Reference_ID")
            cols = ["Virus", "Reference_ID"] + cols
            return combined_df[cols]
        return combined_df

    def _iter_nonempty_files(self) -> Generator[Path | str]:
        """
        Yields input files that exist and are not empty.

        Iterates over the list of input files and yields each file path if the file exists
        and its size is greater than zero bytes.

        Yields
        ------
        Path | str
            Path to an input file that exists and is not empty.
        """
        for infile in self.input_files:
            if os.path.exists(infile) and os.path.getsize(infile) > 0:
                yield infile

    def _combine_coverage(self, file_map: dict[str | Path, tuple[str, str]]) -> None:
        """
        Combine coverage tabular files into a single DataFrame and write to output.

        Reads multiple tabular coverage files, adds virus and reference ID metadata from the provided file map,
        concatenates them into a single DataFrame, and writes the combined result to the specified output file.
        If no non-empty files are found, writes an empty coverage table.

        Parameters
        ----------
        file_map : dict[str | Path, tuple[str, str]]
            A mapping from input file paths to a tuple containing the virus name and reference ID.

        Returns
        -------
        None
            This method writes the combined DataFrame to an output file and does not return a value.
        """
        column_names = [
            "Sample_name",
            "Width_at_mincov_1",
            "Width_at_mincov_5",
            "Width_at_mincov_10",
            "Width_at_mincov_50",
            "Width_at_mincov_100",
        ]
        dfs = []
        for infile in self._iter_nonempty_files():
            df = self._read_tabular_file(infile, names=column_names)
            if df is not None:
                virus, refid = file_map.get(infile, ("Unknown", "Unknown"))
                df["Virus"] = virus
                df["Reference_ID"] = refid
                dfs.append(df)

        if dfs:
            combined_df = pd.concat(dfs, ignore_index=True)
            combined_df = combined_df[
                [
                    "Sample_name",
                    "Virus",
                    "Reference_ID",
                    "Width_at_mincov_1",
                    "Width_at_mincov_5",
                    "Width_at_mincov_10",
                    "Width_at_mincov_50",
                    "Width_at_mincov_100",
                ]
            ]
            combined_df.to_csv(self.output, sep=self.separator, index=False)
        else:
            self._write_empty_tabular("coverage")

    def _combine_mutations(self, file_map: dict[str | Path, tuple[str, str]]) -> None:
        """
        Combine mutation data from multiple tabular files into a single DataFrame and write to output.

        This method iterates over a set of non-empty mutation files, reads each file into a DataFrame,
        annotates each row with the corresponding virus name from the provided file map, and concatenates
        all DataFrames into a single DataFrame. The combined DataFrame is then reordered and written to
        the specified output file. If no non-empty files are found, an empty mutations file is written.

        Parameters
        ----------
        file_map : dict[str | Path, tuple[str, str]]
            A mapping from input file paths to a tuple containing the virus name and an additional string
            (unused in this method). The virus name is added as a column to each DataFrame.

        Returns
        -------
        None
            This method does not return a value. The combined mutations data is written to the output file.
        """
        dfs = []
        for infile in self._iter_nonempty_files():
            df = self._read_tabular_file(infile)
            if df is not None:
                virus, _ = file_map.get(infile, ("Unknown", "Unknown"))
                df["Virus"] = virus
                dfs.append(df)

        if dfs:
            combined_df = pd.concat(dfs, ignore_index=True)
            combined_df = self._reorder_mutations_columns(combined_df)
            combined_df.to_csv(self.output, sep=self.separator, index=False)
        else:
            self._write_empty_tabular("mutations")

    def _combine_amplicon_coverage(self, file_map: dict[str | Path, tuple[str, str]]) -> None:
        """
        Combine amplicon coverage tabular files into a single DataFrame with metadata columns.

        Iterates over non-empty amplicon coverage files, reads each as a DataFrame, and appends
        virus and reference ID metadata columns based on the provided file_map. The resulting
        DataFrames are concatenated, metadata columns are moved to the front, and the combined
        DataFrame is written to the output file. If no valid files are found, an empty CSV is written.

        Parameters
        ----------
        file_map : dict[str | Path, tuple[str, str]]
            A mapping from input file paths to a tuple containing the virus name and reference ID.

        Returns
        -------
        None
            This method writes the combined DataFrame to the output file specified by `self.output`.
        """
        dfs = []
        for infile in self._iter_nonempty_files():
            df = self._read_tabular_file(infile)
            if df is not None:
                virus, refid = file_map.get(infile, ("Unknown", "Unknown"))
                df["Virus"] = virus
                df["Reference_ID"] = refid
                dfs.append(df)

        if dfs:
            combined_df = pd.concat(dfs, ignore_index=True)
            combined_df = self._reorder_amplicon_columns(combined_df)
            combined_df.to_csv(self.output, index=False)
        else:
            pd.DataFrame().to_csv(self.output, index=False)

    def _build_file_map(self) -> dict[str | Path, tuple[str, str]]:
        """
        Builds a mapping from input file paths to tuples containing virus and reference IDs.

        Returns
        -------
        dict[str | Path, tuple[str, str]]
            A dictionary where each key is an input file path (as a string or Path object),
            and each value is a tuple containing the corresponding virus name and reference ID.
        """
        return dict(zip(self.input_files, zip(self.virus_list, self.refid_list)))

    def _combine_tabular(self) -> None:
        """
        Combine tabular files with metadata based on the specified file type.

        Depending on the value of `self.file_type`, this method delegates the combination
        process to the appropriate handler for coverage, mutations, or amplicon coverage
        tabular files. The method first builds a mapping of files to be combined, then
        calls the corresponding combination method.

        Raises
        ------
        ValueError
            If `self.file_type` is not one of the supported types ("coverage", "mutations", "amplicon_coverage").
        """
        file_map = self._build_file_map()
        if self.file_type == "coverage":
            self._combine_coverage(file_map)
        elif self.file_type == "mutations":
            self._combine_mutations(file_map)
        elif self.file_type == "amplicon_coverage":
            self._combine_amplicon_coverage(file_map)


if __name__ == "__main__":
    CombineTabular.main()
