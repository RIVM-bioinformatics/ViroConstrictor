import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PROJECT_ROOT.joinpath("ViroConstrictor/workflow")))
from ViroConstrictor.workflow.main.scripts.extract_gff import ExtractGff  # isort:skip


def test_extract_gff():
    input_gff = PROJECT_ROOT / "tests" / "unit" / "data" / "ESIB_EQA_2024_SARS1_01_features.gff"
    output_gff = PROJECT_ROOT / "tests" / "unit" / "data" / "ESIB_EQA_2024_SARS1_01_features_output.gff"
    ref_id = "NC_045512.2"

    gff_extractor = ExtractGff(input=input_gff, output=output_gff, ref_id=ref_id)
    gff_extractor.run()
