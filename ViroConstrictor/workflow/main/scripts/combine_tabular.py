from argparse import ArgumentParser
from pathlib import Path
from typing import Literal

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
        self._combine_tabular()
    
    def _combine_tabular(self) -> None:
        """Combine tabular files with metadata."""
        file_map = self._build_file_map()
        if self.file_type == "coverage":
            self._combine_coverage(file_map)
        elif self.file_type == "mutations":
            self._combine_mutations(file_map)
        elif self.file_type == "amplicon_coverage":
            self._combine_amplicon_coverage(file_map)

    def _build_file_map(self) -> dict[str | Path, tuple[str, str]]:
        return dict(zip(self.input_files, zip(self.virus_list, self.refid_list)))

    def _iter_nonempty_files(self):
        import os

        for infile in self.input_files:
            if os.path.exists(infile) and os.path.getsize(infile) > 0:
                yield infile

    def _combine_coverage(self, file_map: dict[str | Path, tuple[str, str]]) -> None:
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
            df = self._read_coverage_file(infile, column_names)
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
            self._write_empty_coverage()

    def _read_coverage_file(self, infile, column_names):
        try:
            df = pd.read_csv(infile, sep=self.separator, header=None, names=column_names)
        except pd.errors.EmptyDataError:
            return None
        return df if not df.empty else None

    def _write_empty_coverage(self) -> None:
        header = (
            "Sample_name\tVirus\tReference_ID\tWidth_at_mincov_1\tWidth_at_mincov_5"
            "\tWidth_at_mincov_10\tWidth_at_mincov_50\tWidth_at_mincov_100\n"
        )
        with open(self.output, "w") as f:
            f.write(header)

    def _combine_mutations(self, file_map: dict[str | Path, tuple[str, str]]) -> None:
        dfs = []
        for infile in self._iter_nonempty_files():
            df = self._read_mutations_file(infile)
            if df is not None:
                virus, _ = file_map.get(infile, ("Unknown", "Unknown"))
                df["Virus"] = virus
                dfs.append(df)

        if dfs:
            combined_df = pd.concat(dfs, ignore_index=True)
            combined_df = self._reorder_mutations_columns(combined_df)
            combined_df.to_csv(self.output, sep=self.separator, index=False)
        else:
            self._write_empty_mutations()

    def _read_mutations_file(self, infile):
        try:
            df = pd.read_csv(infile, sep=self.separator)
        except pd.errors.EmptyDataError:
            return None
        return df if not df.empty else None

    def _reorder_mutations_columns(self, combined_df: pd.DataFrame) -> pd.DataFrame:
        cols = combined_df.columns.tolist()
        if "Sample" in cols and "Virus" in cols:
            cols.remove("Virus")
            sample_idx = cols.index("Sample")
            cols.insert(sample_idx + 1, "Virus")
            return combined_df[cols]
        return combined_df

    def _write_empty_mutations(self) -> None:
        header = "Sample\tVirus\tReference_ID\tPosition\tReference_Base\tVariant_Base\tDepth\n"
        with open(self.output, "w") as f:
            f.write(header)

    def _combine_amplicon_coverage(self, file_map: dict[str | Path, tuple[str, str]]) -> None:
        dfs = []
        for infile in self._iter_nonempty_files():
            df = self._read_amplicon_coverage_file(infile)
            if df is not None:
                virus, refid = file_map.get(infile, ("Unknown", "Unknown"))
                df["Virus"] = virus
                df["Reference_ID"] = refid
                dfs.append(df)

        if dfs:
            combined_df = pd.concat(dfs, ignore_index=True)
            combined_df = self._move_amplicon_metadata_front(combined_df)
            combined_df.to_csv(self.output, index=False)
        else:
            pd.DataFrame().to_csv(self.output, index=False)

    def _read_amplicon_coverage_file(self, infile):
        try:
            df = pd.read_csv(infile)
        except pd.errors.EmptyDataError:
            return None
        return None if df.empty else df

    def _move_amplicon_metadata_front(self, combined_df: pd.DataFrame) -> pd.DataFrame:
        cols = combined_df.columns.tolist()
        if "Virus" in cols and "Reference_ID" in cols:
            cols.remove("Virus")
            cols.remove("Reference_ID")
            cols = ["Virus", "Reference_ID"] + cols
            return combined_df[cols]
        return combined_df


if __name__ == "__main__":
    CombineTabular.main()