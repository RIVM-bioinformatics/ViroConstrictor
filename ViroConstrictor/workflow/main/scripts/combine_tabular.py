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
        import os
        
        # Create mapping of file to (virus, refid)
        file_map = dict(zip(self.input_files, zip(self.virus_list, self.refid_list)))
        
        dfs = []
        
        if self.file_type == "coverage":
            # Coverage files have NO HEADER
            column_names = [
                'Sample_name', 'Width_at_mincov_1', 'Width_at_mincov_5',
                'Width_at_mincov_10', 'Width_at_mincov_50', 'Width_at_mincov_100'
            ]
            
            for infile in self.input_files:
                if os.path.exists(infile) and os.path.getsize(infile) > 0:
                    try:
                        df = pd.read_csv(infile, sep=self.separator, header=None, names=column_names)
                        if not df.empty:
                            virus, refid = file_map.get(infile, ('Unknown', 'Unknown'))
                            df['Virus'] = virus
                            df['Reference_ID'] = refid
                            dfs.append(df)
                    except pd.errors.EmptyDataError:
                        continue
            
            if dfs:
                combined_df = pd.concat(dfs, ignore_index=True)
                # Reorder columns
                combined_df = combined_df[[
                    'Sample_name', 'Virus', 'Reference_ID',
                    'Width_at_mincov_1', 'Width_at_mincov_5', 'Width_at_mincov_10',
                    'Width_at_mincov_50', 'Width_at_mincov_100'
                ]]
                combined_df.to_csv(self.output, sep=self.separator, index=False)
            else:
                with open(self.output, 'w') as f:
                    f.write('Sample_name\tVirus\tReference_ID\tWidth_at_mincov_1\tWidth_at_mincov_5\tWidth_at_mincov_10\tWidth_at_mincov_50\tWidth_at_mincov_100\n')
        
        elif self.file_type == "mutations":
            # Mutations files have headers
            for infile in self.input_files:
                if os.path.exists(infile) and os.path.getsize(infile) > 0:
                    try:
                        df = pd.read_csv(infile, sep=self.separator)
                        if not df.empty:
                            virus, refid = file_map.get(infile, ('Unknown', 'Unknown'))
                            df['Virus'] = virus
                            dfs.append(df)
                    except pd.errors.EmptyDataError:
                        continue
            
            if dfs:
                combined_df = pd.concat(dfs, ignore_index=True)
                # Reorder columns to put Virus after Sample
                cols = combined_df.columns.tolist()
                if 'Sample' in cols and 'Virus' in cols:
                    cols.remove('Virus')
                    sample_idx = cols.index('Sample')
                    cols.insert(sample_idx + 1, 'Virus')
                    combined_df = combined_df[cols]
                combined_df.to_csv(self.output, sep=self.separator, index=False)
            else:
                with open(self.output, 'w') as f:
                    f.write('Sample\tVirus\tReference_ID\tPosition\tReference_Base\tVariant_Base\tDepth\n')
        
        elif self.file_type == "amplicon_coverage":
            # Amplicon coverage files are CSVs with headers
            for infile in self.input_files:
                if os.path.exists(infile) and os.path.getsize(infile) > 0:
                    try:
                        df = pd.read_csv(infile)
                        if not df.empty:
                            virus, refid = file_map.get(infile, ('Unknown', 'Unknown'))
                            df['Virus'] = virus
                            df['Reference_ID'] = refid
                            dfs.append(df)
                    except pd.errors.EmptyDataError:
                        continue
            
            if dfs:
                combined_df = pd.concat(dfs, ignore_index=True)
                # Move Virus and Reference_ID to front
                cols = combined_df.columns.tolist()
                if 'Virus' in cols and 'Reference_ID' in cols:
                    cols.remove('Virus')
                    cols.remove('Reference_ID')
                    cols = ['Virus', 'Reference_ID'] + cols
                    combined_df = combined_df[cols]
                combined_df.to_csv(self.output, index=False)
            else:
                pd.DataFrame().to_csv(self.output, index=False)


if __name__ == "__main__":
    CombineTabular.main()