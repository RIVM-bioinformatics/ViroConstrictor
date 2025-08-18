from pathlib import Path
from typing import Any, Protocol, Type

import pytest

from ViroConstrictor.workflow.scripts import (
    AmpliconCovs,
    Boc,
    Clipper,
    ConcatAmpliconCovs,
    CountMappedReads,
    ExtractGff,
    FilterBed,
    FilterBestMatchingRef,
    FilterGff,
    FilterReferences,
    GroupAminoAcids,
    GroupRefs,
    PrepareRefs,
    VcfToTsv,
)


class ScriptProtocol(Protocol):
    def __init__(self, **kwargs: Any) -> None: ...

    def run(self) -> None: ...


SCRIPTS: list[dict[str, Any]] = [
    {
        "class": Boc,
        "args": {
            "input": "tests/unit/data/input_boc.fasta",
            "output": "tests/unit/data/output_boc.txt",
            "samplename": "sample1",
            "coverages": "tests/unit/data/coverage.tsv",
        },
    },
    {
        "class": ExtractGff,
        "args": {
            "input": "tests/unit/data/input_gff.gff",
            "output": "tests/unit/data/output_gff.gff",
            "refID": "ref1",
        },
    },
    {
        "class": AmpliconCovs,
        "args": {
            "input": "tests/unit/data/input_amplicon_covs.fasta",
            "output": "tests/unit/data/output_amplicon_covs.txt",
        },
    },
    {
        "class": Clipper,
        "args": {
            "input": "tests/unit/data/input_clipper.fasta",
            "output": "tests/unit/data/output_clipper.txt",
        },
    },
    {
        "class": ConcatAmpliconCovs,
        "args": {
            "input": ["tests/unit/data/input_amplicon_cov1.txt", "tests/unit/data/input_amplicon_cov2.txt"],
            "output": "tests/unit/data/output_concat_amplicon_covs.txt",
        },
    },
    {
        "class": FilterBed,
        "args": {
            "input": "tests/unit/data/input_bed.bed",
            "refdata": "tests/unit/data/ref_data.csv",
            "output": "tests/unit/data/output_bed.bed",
            "updatedstats": "tests/unit/data/updated_stats.csv",
        },
    },
    {
        "class": FilterBestMatchingRef,
        "args": {
            "input": "tests/unit/data/input_counts.csv",
            "inputref": "tests/unit/data/input_refs.fasta",
            "filtref": "tests/unit/data/output_filtered_refs.fasta",
            "output": "tests/unit/data/output_filtered_counts.csv",
        },
    },
    {
        "class": FilterGff,
        "args": {
            "input": "tests/unit/data/input_gff.gff",
            "refdata": "tests/unit/data/ref_data.csv",
            "output": "tests/unit/data/output_filtered_gff.gff",
            "updatedstats": "tests/unit/data/updated_stats.csv",
        },
    },
    {
        "class": GroupAminoAcids,
        "args": {
            "input": "tests/unit/data/input_amino_acids.fasta",
            "output": "tests/unit/data/output_grouped_amino_acids.fasta",
        },
    },
    {
        "class": GroupRefs,
        "args": {
            "input_refs": ["tests/unit/data/input_ref1.fasta", "tests/unit/data/input_ref2.fasta"],
            "input_stats": ["tests/unit/data/input_stats1.csv", "tests/unit/data/input_stats2.csv"],
            "output_ref": "tests/unit/data/output_grouped_refs.fasta",
            "output_stats": "tests/unit/data/output_grouped_stats.csv",
            "sample": "sample1",
        },
    },
    {
        "class": PrepareRefs,
        "args": {
            "input": "tests/unit/data/input_prepare_refs.fasta",
            "output": "tests/unit/data/output_prepare_refs.fasta",
        },
    },
    {
        "class": VcfToTsv,
        "args": {
            "input": "tests/unit/data/input_vcf.vcf",
            "output": "tests/unit/data/output_vcf_to_tsv.tsv",
        },
    },
]


MATCH_REF_SCRIPTS: list[dict[str, Any]] = [
    {
        "class": FilterReferences,
        "args": {
            "input": "tests/unit/data/input_refs.fasta",
            "output": "tests/unit/data/output_refs.fasta",
            "wildcard_segment": "segment1",
        },
    },
    {
        "class": CountMappedReads,
        "args": {
            "input": "tests/unit/data/input_bam.bam",
            "output": "tests/unit/data/output_counts.csv",
        },
    },
    {
        "class": FilterBed,
        "args": {
            "input": "tests/unit/data/input_bed.bed",
            "refdata": "tests/unit/data/ref_data.csv",
            "output": "tests/unit/data/output_bed.bed",
            "updatedstats": "tests/unit/data/updated_stats.csv",
        },
    },
    {
        "class": FilterBestMatchingRef,
        "args": {
            "input": "tests/unit/data/input_counts.csv",
            "inputref": "tests/unit/data/input_refs.fasta",
            "filtref": "tests/unit/data/output_filtered_refs.fasta",
            "output": "tests/unit/data/output_filtered_counts.csv",
        },
    },
    {
        "class": FilterGff,
        "args": {
            "input": "tests/unit/data/input_gff.gff",
            "refdata": "tests/unit/data/ref_data.csv",
            "output": "tests/unit/data/output_filtered_gff.gff",
            "updatedstats": "tests/unit/data/updated_stats.csv",
        },
    },
    {
        "class": GroupRefs,
        "args": {
            "input_refs": ["tests/unit/data/input_ref1.fasta", "tests/unit/data/input_ref2.fasta"],
            "input_stats": ["tests/unit/data/input_stats1.csv", "tests/unit/data/input_stats2.csv"],
            "output_ref": "tests/unit/data/output_grouped_refs.fasta",
            "output_stats": "tests/unit/data/output_grouped_stats.csv",
            "sample": "sample1",
        },
    },
]

all_scripts = SCRIPTS + MATCH_REF_SCRIPTS
all_script_ids = [f"{script['class__.__name__']}_{i}" for i, script in enumerate(all_scripts)]


@pytest.mark.parametrize("script_data", all_scripts, ids=all_script_ids)
def test_scripts(script_data: dict[str, Any], tmp_path: Path) -> None:
    """
    Generalized test for all scripts.

    Parameters
    ----------
    script_data : dict[str, Any]
        A dictionary containing the script class and its arguments.
    tmp_path : Path
        Temporary directory provided by pytest for test outputs.
    """
    # Prepare arguments, replacing output paths with temporary paths
    args: dict[str, Any] = script_data["args"].copy()
    for key, value in args.items():
        if key == "output" or key.endswith("output"):
            args[key] = tmp_path / Path(value).name

    # Instantiate and run the script
    script_class: Type[ScriptProtocol] = script_data["class"]
    script_instance = script_class(**args)
    script_instance.run()

    # Assert that the output file was created
    output_path: Path = args["output"]
    assert output_path.exists(), f"Output file {output_path} was not created."
