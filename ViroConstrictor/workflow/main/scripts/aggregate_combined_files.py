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
        self._aggregate_files()
    
    def _aggregate_files(self) -> None:
        """Concatenate already-processed files with headers."""
        import os
        
        dfs = []
        
        # Read all input files
        for infile in self.input_files:
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
        
        # Write combined output
        if dfs:
            combined_df = pd.concat(dfs, ignore_index=True)
            if self.file_type == "amplicon_coverage":
                combined_df.to_csv(self.output, index=False)
            else:
                combined_df.to_csv(self.output, sep=self.separator, index=False)
        else:
            # Write empty file with appropriate header
            self._write_empty_file()
    
    def _write_empty_file(self) -> None:
        """Write an empty file with appropriate headers based on file type."""
        headers = {
            "mutations": "Sample\tVirus\tReference_ID\tPosition\tReference_Base\tVariant_Base\tDepth\n",
            "coverage": "Sample_name\tVirus\tReference_ID\tWidth_at_mincov_1\tWidth_at_mincov_5\tWidth_at_mincov_10\tWidth_at_mincov_50\tWidth_at_mincov_100\n",
            "amplicon_coverage": "",  # Empty CSV
        }
        
        with open(self.output, 'w') as f:
            if self.file_type in headers:
                f.write(headers[self.file_type])


if __name__ == "__main__":
    AggregateCombinedFiles.main()
