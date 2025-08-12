import difflib
import inspect
import json
import os
import re
from collections import defaultdict
from typing import Any, List, Tuple

from rich import print

from ViroConstrictor.logging import log

# Load the preset aliases from a JSON file (preset_aliases.json)
aliases = json.load(
    open(
        os.path.join(os.path.abspath(os.path.dirname(__file__)), "preset_aliases.json")
    )
)

# Load the preset parameters from a JSON file (preset_params.json)
presets = json.load(
    open(os.path.join(os.path.abspath(os.path.dirname(__file__)), "preset_params.json"))
)


def get_key_from_value(d: dict, value: str) -> str | None:
    """This function finds the key in a dictionary which has a value matching the input value.

    Parameters
    ----------
    d
        The dictionary to search.
    value
        The value to search for.

    Returns
    -------
    The key from the dictionary which has the input value.

    """
    for k, v in d.items():
        if value in v:
            return k


def match_preset_name(targetname: str, use_presets: bool) -> Tuple[str, float]:
    """The function takes a target name and a boolean flag as input, and returns a tuple containing the
    best matching preset name and a score based on string similarity, or a default value if the flag is
    False or no match is found with a high enough similarity score.

    Parameters
    ----------
    targetname : str
        The name of the target that needs to be matched with a preset name.
    use_presets : bool
        A boolean parameter that determines whether to use preset values or not. If it is set to True, the
    function will use preset values to match the targetname. If it is set to False, the function will
    return a default value.

    Returns
    -------
        A tuple containing a string and a float. The string represents the matched preset name, and the
    float represents the similarity score between the target name and the matched preset name. If the
    use_presets parameter is False, the function returns the string "DEFAULT" and the float 0. If no
    match is found with a similarity score greater than or equal to 0.4, the function returns the string
    "DEFAULT" and the float 0. If a match is found with a similarity score greater than or equal to 0.4,
    the function returns the matched preset name and the similarity score.
    """
    if not use_presets:
        return "DEFAULT", float(0)
    # regex to remove all special characters from targetname except underscores and dashes
    query = re.sub(r"[^_a-zA-Z0-9/-]+", "", targetname).upper()

    if query == "DEFAULT":
        return "DEFAULT", float(1)

    # flatten list of lists aliases.values() into a single list
    aliases_list = [item for sublist in aliases.values() for item in sublist]

    best_match = difflib.get_close_matches(query, aliases_list, cutoff=0.0, n=1)[0]
    score = difflib.SequenceMatcher(None, a=query, b=best_match).ratio()

    if score < 0.40:
        return "DEFAULT", float(0)
    if matched_preset := get_key_from_value(aliases, best_match):
        return matched_preset, score
    return "DEFAULT", float(0)


def collapse_preset_group(
    preset_name: str, stages: List[str], stage_identifier: str
) -> dict[str, str]:
    """
    Collapse one or more preset groups into a single dictionary.

    Parameters
    ----------
    preset_name : str
        The name of the preset group.
    stages : list of str
        The list of stages, or subgroups, to include in the collapsed dictionary.
    stage_identifier : str
        The identifier to prepend to the keys in the collapsed dictionary.

    Returns
    -------
    collapsed_dict : dict[str, str]
        The collapsed dictionary.

    Examples
    --------
    >>> presets = {
    ...     "group1": {
    ...         "stage1": {"key1": "value1", "key2": "value2"},
    ...         "stage2": {"key3": "value3", "key4": "value4"}
    ...     },
    ...     "group2": {
    ...         "stage1": {"key5": "value5", "key6": "value6"},
    ...         "stage2": {"key7": "value7", "key8": "value8"}
    ...     }
    ... }
    >>> collapse_preset_group("group1", ["stage1", "stage2"], "prefix")
    {'prefix_key1': 'value1', 'prefix_key2': 'value2', 'prefix_key3': 'value3', 'prefix_key4': 'value4'}
    """
    temp_dict: dict[str, str] = {}
    for k, v in presets[preset_name].items():
        if k in stages:
            for nk, nv in v.items():
                temp_dict[f"{stage_identifier}_{nk}"] = nv
    return temp_dict


def get_preset_parameter(preset_name: str, parameter_name: str) -> str:
    """
    Flexibly get predefined tool-parameters from one or more presets.

    Parameters
    ----------
    preset_name : str
        The name of the preset.
    parameter_name : str
        The name of the parameter.

    Returns
    -------
    str
        The value of the parameter.

    Raises
    ------
    KeyError
        If the preset or parameter is not found.

    Notes
    -----
    This function retrieves the value of a parameter from a preset. The preset can have different values for different stages of execution. The stages are identified by a stage identifier, which is obtained dynamically using the `inspect` module. The stage identifier is stored in the global variable `VC_STAGE`.

    The function first collapses the preset groups into dictionaries for the main and matchref stages. This results in dictionaries with only the override values for the specific stage that is being called.

    The dictionaries for the main and matchref stages are then merged together, with the main stage taking precedence over the matchref stage. This results in one dictionary with the override values, where the stage identifier is prepended to the key.

    The function uses a `defaultdict` to allow for a default value of `None` if the key is not found in the dictionary. This indicates that there is no override value and the default value should be used.

    If the parameter is a dictionary, it means that the parameter is not found in the specific preset group. In this case, the function fetches the parameter from the default preset. The default preset is obtained by collapsing the preset group for the "DEFAULT" stage.

    Examples
    --------
    >>> get_preset_parameter("preset1", "parameter1")
    'value1'
    >>> get_preset_parameter("preset2", "parameter2")
    'value2'

    """

    # Using an inspect here to have a dynamic stage identifier that cannot be changed.
    # Additionally, this allows for the function to be used in any stage without needing to pass the stage identifier at every function call.
    _stage_identifier = inspect.stack()[1][0].f_globals["VC_STAGE"]

    # collapse the preset groups into dictionaries for the main and matchref stages
    # This results in the dictionaries with only the override values for the specific stage that is being called.
    presets_main = collapse_preset_group(
        preset_name, ["STAGE_MAIN", "STAGE_GLOBAL"], _stage_identifier
    )
    presets_matchref = collapse_preset_group(
        preset_name, ["STAGE_MATCHREF", "STAGE_GLOBAL"], _stage_identifier
    )

    # merge the two dictionaries together, with the main stage taking precedence over the matchref stage
    # this results in one dictionary with the override values with the stage prepended to the key
    # We're using a defaultdict here to allow for the default value to be None if the key is not found, thus indicating that there is no override value and the default value should be used.
    preset = defaultdict(dict[str, str], presets_matchref | presets_main, default=None)

    # fetch the desired parameter, returns either the override value or None if the key does not exist.
    parameter: Any = preset[f"{_stage_identifier}_{parameter_name}"]

    # if the parameter is a dictionary, it means that the parameter is not found in the specific preset group
    if isinstance(parameter, dict):
        # fetch the parameter from the default preset
        preset_params = collapse_preset_group(
            "DEFAULT", ["STAGE_MAIN", "STAGE_GLOBAL"], _stage_identifier
        )
        log.debug(
            f"Python code :: Preset :: {preset_name} specific parameter '[yellow]{_stage_identifier}_{parameter_name}[/yellow]' not found, inheriting parameter value from 'DEFAULT'."
        )

        return preset_params[f"{_stage_identifier}_{parameter_name}"]

    log.debug(
        f"Python code :: Preset :: Using {preset_name} specific parameter '[yellow]{_stage_identifier}_{parameter_name}[/yellow]'"
    )
    return parameter
