import difflib
import importlib.resources
import inspect
import json
import re
from collections import defaultdict
from functools import cache
from typing import Any, List, Tuple

from rich import print

from ViroConstrictor.logging import log


def _load_preset_resource(resource_name: str) -> dict:
    """Load a single preset JSON resource from viroconstrictor-data package.

    Parameters
    ----------
    resource_name : str
        JSON file name inside the ``viroconstrictor_data/presets`` resource directory.

    Returns
    -------
    dict
        Parsed JSON content for the requested preset resource.

    Raises
    ------
    ImportError
        If the viroconstrictor-data package is not installed, or if required resources are missing/invalid.

    """
    try:
        pkg = importlib.resources.files("viroconstrictor_data") / "presets"
        return json.loads((pkg / resource_name).read_text(encoding="utf-8"))
    except (ModuleNotFoundError, json.JSONDecodeError, OSError):
        raise ImportError("viroconstrictor-data preset resources are missing or invalid. Run: pip install --upgrade --force-reinstall viroconstrictor-data") from None


# This function is currently not used but it's there for future reference if we want to load both presets and aliases at the same time.
def _load_preset_data() -> tuple[dict, dict]:
    """Load preset parameter and alias dictionaries from viroconstrictor-data.

    Returns
    -------
    tuple[dict, dict]
        A two-item tuple containing preset parameters and preset aliases.

    """
    return _load_preset_resource("preset_params.json"), _load_preset_resource("preset_aliases.json")


@cache
def _get_presets() -> dict:
    """Load and cache preset parameters.

    Returns
    -------
    dict
        Preset parameter definitions loaded from ``preset_params.json``.

    """
    return _load_preset_resource("preset_params.json")


@cache
def _get_aliases() -> dict:
    """Load and cache preset aliases.

    Returns
    -------
    dict
        Preset alias mappings loaded from ``preset_aliases.json``.

    """
    return _load_preset_resource("preset_aliases.json")


def get_key_from_value(d: dict, value: str) -> str | None:
    """Find the key in a dictionary containing a matching value.

    Parameters
    ----------
    d : dict
        The dictionary to search.
    value : str
        The value to search for.

    Returns
    -------
    str or None
        The key from the dictionary with the input value, or None if not found.

    Examples
    --------
    Find a key containing a specific value in a dictionary:

    >>> data = {"virus": ["SARSCOV2", "COVID-19"], "platform": ["Nanopore", "Illumina"]}
    >>> get_key_from_value(data, "SARSCOV2")
    'virus'

    Return None when the value is not found:

    >>> get_key_from_value(data, "NonExistent") is None
    True

    """
    for k, v in d.items():
        if value in v:
            return k


def match_preset_name(targetname: str, use_presets: bool) -> Tuple[str, float]:
    """Match a target name to the best matching preset using fuzzy string matching.

    Parameters
    ----------
    targetname : str
        The name of the target that needs to be matched with a preset name.
    use_presets : bool
        If True, match the targetname to available presets. If False, return DEFAULT with score 0.

    Returns
    -------
    tuple[str, float]
        A tuple containing the preset name and similarity score. If no match with similarity ≥ 0.4 is
        found or use_presets is False, returns ("DEFAULT", 0.0).

    Examples
    --------
    Match a target name with fuzzy matching enabled:

    >>> preset, score = match_preset_name("sarscov2", use_presets=True)
    >>> preset
    'SARSCOV2'
    >>> score
    1.0

    Fuzzy match with typo or variant spelling:

    >>> preset, score = match_preset_name("SARS-COV-2", use_presets=True)
    >>> preset
    'SARSCOV2'
    >>> score  # doctest: +SKIP
    0.95  # High similarity but not perfect match

    Return DEFAULT when presets are disabled:

    >>> preset, score = match_preset_name("SARSCOV2", use_presets=False)
    >>> (preset, score)
    ('DEFAULT', 0.0)

    """
    if not use_presets:
        return "DEFAULT", float(0)
    # regex to remove all special characters from targetname except underscores and dashes
    query = re.sub(r"[^_a-zA-Z0-9/-]+", "", targetname).upper()

    if query == "DEFAULT":
        return "DEFAULT", float(1)

    preset_aliases = _get_aliases()

    # flatten list of lists aliases.values() into a single list
    aliases_list = [item for sublist in preset_aliases.values() for item in sublist]

    best_match = difflib.get_close_matches(query, aliases_list, cutoff=0.0, n=1)[0]
    score = difflib.SequenceMatcher(None, a=query, b=best_match).ratio()

    if score < 0.40:
        return "DEFAULT", float(0)
    if matched_preset := get_key_from_value(preset_aliases, best_match):
        return matched_preset, score
    return "DEFAULT", float(0)


def collapse_preset_group(preset_name: str, stages: List[str], stage_identifier: str) -> dict[str, str]:
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
    preset_data = _get_presets()
    temp_dict: dict[str, str] = {}
    for k, v in preset_data[preset_name].items():
        if k in stages:
            for nk, nv in v.items():
                temp_dict[f"{stage_identifier}_{nk}"] = nv
    return temp_dict


def get_preset_parameter(preset_name: str, parameter_name: str, stage_identifier: str = "") -> Any:
    """Retrieve predefined tool parameters from a preset, with fallback to DEFAULT preset.

    Dynamically retrieves preset parameters based on the current workflow stage (MAIN or MATCHREF).
    The stage identifier is obtained from the calling function's global `VC_STAGE` variable unless
    explicitly provided.

    Parameters
    ----------
    preset_name : str
        The name of the preset to retrieve parameters from.
    parameter_name : str
        The name of the parameter to retrieve.
    stage_identifier : str, optional
        Explicit stage identifier (default: auto-detect from caller's `VC_STAGE`).

    Returns
    -------
    Any
        The parameter value from the preset. Falls back to the DEFAULT preset if the parameter is
        not found in the specified preset.

    Raises
    ------
    KeyError
        If the parameter is not found in either the specified preset or DEFAULT preset.

    Notes
    -----
    **Stage Identifier Detection**
        If stage_identifier is not provided, the function dynamically obtains it from the caller's
        global `VC_STAGE` variable using the `inspect` module. This ensures the function can be used
        in any stage without requiring explicit parameter passing.

    **Preset Resolution Logic**
        1. Collapses preset groups into dictionaries for MAIN and MATCHREF stages
        2. Merges dictionaries with MAIN stage taking precedence over MATCHREF
        3. Uses a defaultdict to return None for missing keys, indicating no override exists
        4. Prepends stage identifier to parameter keys (e.g., "MAIN_AlignmentPreset")

    **Fallback Mechanism**
        If a parameter is not found in the specified preset (indicated by a dict-type value),
        the function automatically falls back to the DEFAULT preset. This ensures parameters
        are always available and prevents KeyError exceptions for optional overrides.

    **Key Format**
        Parameters are stored with stage identifier prefix:
        - STAGE_MAIN: Parameters for main workflow stage
        - STAGE_MATCHREF: Parameters for match-reference stage
        - STAGE_GLOBAL: Parameters shared across all stages

    Examples
    --------
    Retrieve Minimap2 base alignment settings for a SARSCOV2 preset in the MAIN stage:

    >>> settings = get_preset_parameter("SARSCOV2", "Minimap2_Settings_Base")
    >>> # Returns: "--secondary=no", retrieved from DEFAULT - GLOBAL since SARSCOV2 does not override this parameter

    Retrieve adapter removal preset specific to platform:

    >>> adapters = get_preset_parameter("INFLUENZA", "FastP_AdapterRemoval_Settings_illumina")
    >>> # Automatically falls back to DEFAULT preset if INFLUENZA-specific setting not found

    Use explicit stage identifier for testing or non-standard workflows:

    >>> params = get_preset_parameter("SARSCOV2", "Clipper_FilterParams_nanopore", stage_identifier="MAIN")
    >>> # Explicitly specifies MAIN stage instead of auto-detecting from VC_STAGE

    """

    # Using an inspect here to have a dynamic stage identifier that cannot be changed.
    # Additionally, this allows for the function to be used in any stage without needing to pass the stage identifier at every function call.
    if not stage_identifier:
        _stage_identifier = inspect.stack()[1][0].f_globals["VC_STAGE"]
    else:
        _stage_identifier = stage_identifier

    # collapse the preset groups into dictionaries for the main and matchref stages
    # This results in the dictionaries with only the override values for the specific stage that is being called.
    presets_main = collapse_preset_group(preset_name, ["STAGE_MAIN", "STAGE_GLOBAL"], _stage_identifier)
    presets_matchref = collapse_preset_group(preset_name, ["STAGE_MATCHREF", "STAGE_GLOBAL"], _stage_identifier)

    # merge the two dictionaries together, with the main stage taking precedence over the matchref stage
    # this results in one dictionary with the override values with the stage prepended to the key
    # We're using a defaultdict here to allow for the default value to be None if the key is not found, thus indicating that there is no override value and the default value should be used.
    preset = defaultdict(dict[str, str], presets_matchref | presets_main, default=None)

    # fetch the desired parameter, returns either the override value or None if the key does not exist.
    parameter: Any = preset[f"{_stage_identifier}_{parameter_name}"]

    # if the parameter is a dictionary, it means that the parameter is not found in the specific preset group
    if isinstance(parameter, dict):
        # fetch the parameter from the default preset
        preset_params = collapse_preset_group("DEFAULT", ["STAGE_MAIN", "STAGE_GLOBAL"], _stage_identifier)
        log.debug(
            f"Parameter handling :: Preset :: {preset_name} specific parameter '[yellow]{_stage_identifier}_{parameter_name}[/yellow]' not found, inheriting parameter value from 'DEFAULT'."
        )

        return preset_params[f"{_stage_identifier}_{parameter_name}"]

    log.debug(f"Parameter handling :: Preset :: Using {preset_name} specific parameter '[yellow]{_stage_identifier}_{parameter_name}[/yellow]'")
    return parameter
