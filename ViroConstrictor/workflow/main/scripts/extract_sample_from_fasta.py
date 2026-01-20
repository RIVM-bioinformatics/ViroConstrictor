from argparse import ArgumentParser
from pathlib import Path

from Bio import SeqIO
from helpers.base_script_class import BaseScript


class ExtractSampleFromFasta(BaseScript):
    """Extract sequences for a specific sample from FASTA files."""
    
    def __init__(
        self,
        input: str,
        output: Path | str,
        input_files: list[Path | str],
        sample_name: str,
    ) -> None:
        super().__init__(input, output)
        self.input_files = input_files
        self.sample_name = sample_name
    
    @classmethod
    def add_arguments(cls, parser: ArgumentParser) -> None:
        super().add_arguments(parser)
        parser.add_argument(
            "--input_files",
            nargs="+",
            help="List of input FASTA files to extract from.",
            type=str,
            required=True,
        )
        parser.add_argument(
            "--sample_name",
            help="Sample name to extract sequences for.",
            type=str,
            required=True,
        )
    
    def run(self) -> None:
        self._extract_sequences()
    
    def _extract_sequences(self) -> None:
        """Extract sequences that belong to the specified sample."""
        import os
        
        extracted_records = []
        
        for infile in self.input_files:
            if not os.path.exists(infile) or os.path.getsize(infile) == 0:
                continue
            
            for record in SeqIO.parse(infile, "fasta"):
                # Check if the sequence ID starts with the sample name
                # Format is: sampleID.featureName
                seq_id_parts = record.id.split(".")
                if seq_id_parts[0] == self.sample_name:
                    extracted_records.append(record)
        
        # Write extracted records to output
        with open(self.output, 'w') as out_handle:
            SeqIO.write(extracted_records, out_handle, "fasta")


if __name__ == "__main__":
    ExtractSampleFromFasta.main()