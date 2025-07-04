# pylint: disable=C0103

"""
Construct and write configuration files for SnakeMake
"""

import os

import yaml


# TODO: find another place for this function
def WriteYaml(data: dict, filepath: str) -> str:
    """WriteYaml takes a dictionary and a filepath, and writes the given dictionary to the filepath as a yaml
    file

    Parameters
    ----------
    data : dict
        The data to be written to the file.
    filepath : str
        The path to the file you want to write to.

    Returns
    -------
    filepath : str
        The filepath

    """
    if not os.path.exists(os.path.dirname(filepath)):
        os.makedirs(os.path.dirname(filepath))
    with open(filepath, "w") as file:
        yaml.dump(data, file, default_flow_style=False)
    return filepath
