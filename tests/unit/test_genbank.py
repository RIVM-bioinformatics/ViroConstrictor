from pathlib import Path

from ViroConstrictor.workflow.scripts.genbank import main

PATH = Path("tests/unit/data/test_reference_double.gb")


def test_split_genbank():
    args = [
        str(PATH),
        "--include-target",
    ]
    main(args)
