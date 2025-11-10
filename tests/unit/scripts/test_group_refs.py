import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PROJECT_ROOT.joinpath("ViroConstrictor/workflow")))
from ViroConstrictor.workflow.match_ref.scripts.group_refs import GroupRefs  # isort:skip


def test_group_refs():
    _input = "empty"
    _input_refs = [
        PROJECT_ROOT / "tests" / "unit" / "data" / "group_refs" / "HA_INF1024_B34_4312201729_squiggle_best_ref.fasta",
        PROJECT_ROOT / "tests" / "unit" / "data" / "group_refs" / "MP_INF1024_B34_4312201729_squiggle_best_ref.fasta",
    ]
    _input_stats = [
        PROJECT_ROOT / "tests" / "unit" / "data" / "group_refs" / "HA_INF1024_B34_4312201729_squiggle_best_ref.csv",
        PROJECT_ROOT / "tests" / "unit" / "data" / "group_refs" / "MP_INF1024_B34_4312201729_squiggle_best_ref.csv",
    ]

    output = PROJECT_ROOT / "tests" / "unit" / "data" / "group_refs" / "output" / "output_refs.fasta"
    output_stats = PROJECT_ROOT / "tests" / "unit" / "data" / "group_refs" / "output" / "output_stats.csv"
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
