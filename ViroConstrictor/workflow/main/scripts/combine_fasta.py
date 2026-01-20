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
        self._combine_fasta()
    
    def _combine_fasta(self) -> None:
        """Combine FASTA files with modified headers."""
        import os
        
        # Create mapping of file to (virus, refid)
        file_map = dict(zip(self.input_files, zip(self.virus_list, self.refid_list)))
        
        with open(self.output, 'w') as out_handle:
            for infile in self.input_files:
                if not os.path.exists(infile) or os.path.getsize(infile) == 0:
                    continue
                
                virus, refid = file_map[infile]
                
                for record in SeqIO.parse(infile, "fasta"):
                    # Original header format: >sampleID mincov=X
                    # New format: >sampleID Virus RefID mincov=X
                    header_parts = record.description.split()
                    if len(header_parts) >= 2 and 'mincov=' in header_parts[-1]:
                        sample_id = header_parts[0]
                        mincov_part = header_parts[-1]
                        new_header = f"{sample_id} {virus} {refid} {mincov_part}"
                        record.id = sample_id
                        record.description = new_header
                    else:
                        # Fallback: just append Virus and RefID
                        new_header = f"{record.description} {virus} {refid}"
                        record.description = new_header
                    
                    SeqIO.write(record, out_handle, "fasta")


if __name__ == "__main__":
    CombineFasta.main()