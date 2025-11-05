import sys
from pathlib import Path

project_root = Path(__file__).parent.parent.parent.parent
sys.path.insert(0, str(project_root.joinpath("ViroConstrictor/workflow")))
from ViroConstrictor.workflow.main.scripts.extract_gff import ExtractGff  # isort:skip


def test_extract_gff():
    input_gff = "tests/unit/data/ESIB_EQA_2024_SARS1_01_features.gff"
    output_gff = "tests/unit/data/ESIB_EQA_2024_SARS1_01_features_output.gff"
    ref_id = "NC_045512.2"

    gff_extractor = ExtractGff(input=input_gff, output=output_gff, ref_id=ref_id)
    gff_extractor.run()
