from 

SCRIPTS = [
    {
        "class": GroupAminoAcids,
        "args": {
            "input": "tests/unit/data/aa.faa",
            "output": "tests/unit/data/aa_output.faa",
            "space": "tests/unit/data/sampleinfo.pkl",
        },
    },
    {
        "class": ExtractGff,
        "args": {
            "input": "tests/unit/data/example.gff",
            "output": "tests/unit/data/example_output.gff",
            "ref_id": "example_ref",
        },
    },
    # Add more scripts here
]


@pytest.mark.parametrize("script_data", SCRIPTS)
def test_scripts(script_data, tmp_path: Path) -> None:
    """
    Generalized test for all scripts.

    Parameters
    ----------
    script_data : dict
        A dictionary containing the script class and its arguments.
    tmp_path : Path
        Temporary directory provided by pytest for test outputs.
    """
    # Prepare arguments, replacing output paths with temporary paths
    args = script_data["args"].copy()
    args["output"] = tmp_path / Path(args["output"]).name

    # Instantiate and run the script
    script_instance = script_data["class"](**args)
    script_instance.run()

    # Assert that the output file was created
    assert args["output"].exists(), f"Output file {args['output']} was not created."
