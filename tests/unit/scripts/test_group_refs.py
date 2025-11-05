import sys
from pathlib import Path

project_root = Path(__file__).parent.parent.parent.parent
sys.path.insert(0, str(project_root.joinpath("ViroConstrictor/workflow")))
from ViroConstrictor.workflow.match_ref.scripts.group_refs import GroupRefs  # isort:skip


def test_group_refs():
    _input = "empty"
    _input_refs = [
        "tests/unit/data/group_refs/HA_INF1024_B34_4312201729_squiggle_best_ref.fasta",
        "tests/unit/data/group_refs/MP_INF1024_B34_4312201729_squiggle_best_ref.fasta",
    ]
    _input_stats = [
        "tests/unit/data/group_refs/HA_INF1024_B34_4312201729_squiggle_best_ref.csv",
        "tests/unit/data/group_refs/MP_INF1024_B34_4312201729_squiggle_best_ref.csv",
    ]

    output = "tests/unit/data/group_refs/output/output_refs.fasta"
    output_stats = "tests/unit/data/group_refs/output/output_stats.csv"
    sample_name = "INF1024_B34_4312201729_squiggle"

    grouper = GroupRefs(
        input=_input,
        input_refs=_input_refs,
        input_stats=_input_stats,
        output=output,
        output_stats=output_stats,
        sample=sample_name,
    )
    grouper.run()
