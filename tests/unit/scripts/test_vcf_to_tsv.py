from pathlib import Path

import pandas as pd

from ViroConstrictor.workflow.main.scripts.vcf_to_tsv import VcfToTsv


def test_vcf_to_tsv(tmp_path: Path) -> None:
    input_path = tmp_path / "example.vcf"
    output_path = tmp_path / "result.tsv"
    samplename = "Sample1"

    # Create a mock VCF file
    input_path.write_text(
        """##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t12345\t.\tA\tT\t.\tPASS\tDP=100
chr1\t67890\t.\tG\tN\t.\tPASS\tDP=50
"""
    )

    script = VcfToTsv(input_path, output_path, samplename)
    script.run()

    # Verify the output
    df = pd.read_csv(output_path, sep="\t", header=None)
    assert df.iloc[0, 0] == "Sample1"
    assert df.iloc[0, 1] == "chr1"
    assert df.iloc[0, 2] == 12345
    assert df.iloc[0, 3] == "A"
    assert df.iloc[0, 4] == "T"
    assert df.iloc[0, 5] == 100
