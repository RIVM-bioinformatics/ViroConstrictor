from pathlib import Path

from ViroConstrictor.workflow.scripts.genbank import split_genbank

PATH = Path("tests/unit/data/test_reference_double.gb")


def test_split_genbank():
    split_genbank(PATH, include_target=True)
