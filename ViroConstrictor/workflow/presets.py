import difflib
import re
from typing import Tuple

aliases = {
    "SARSCOV2": [
        "SARS-COV-2",
        "SARS2",
        "COVID",
        "COV",
        "SARSCOV",
        "SARSCOV2",
        "CORONAVIRUS",
    ],
    "INFLUENZA": [
        "FLU",
        "INFLU",
        "INFLUENZA",
        "INFLUENZA_A",
        "INFLUENZA_B",
        "INF",
    ],
    # "MEASLES": ["MEASLES", "MEV", "MEASLES_VIRUS"],
    # "HPV": ["PAPILLOMA_VIRUS", "HPV"]
}

presets = {
    "DEFAULT": {
        "MatchRef_AdditionalAlignmentSettings": "--secondary=no",
        "MatchRef_AdditionalBamFilters": "-q 30",
        "RawAlign_AdditionalSettings": "",
        "BaseBAMFilters": "-F 256 -F 512 -F 4 -F 2048",
        "ClipperSettings_nanopore": "",
        "ClipperSettings_illumina": "",
        "ClipperSettings_iontorrent": "",
        "Fastp_PhredScore_cutoff_illumina": 20,
        "Fastp_PhredScore_cutoff_nanopore": 7,
        "Fastp_PhredScore_cutoff_iontorrent": 20,
        "Fastp_WindowSize_illumina": 5,
        "Fastp_WindowSize_nanopore": 20,
        "Fastp_WindowSize_iontorrent": 15,
        "Fastp_MinReadLength": 100,
        "CleanAlign_AdditionalSettings_nanopore": "",
        "CleanAlign_AdditionalSettings_illumina": "",
        "CleanAlign_AdditionalSettings_iontorrent": "",
    },
    "SARSCOV2": {
        "MatchRef_AdditionalAlignmentSettings": "--secondary=no",
        "MatchRef_AdditionalBamFilters": "-q 50",
        "RawAlign_AdditionalSettings": "",
        "BaseBAMFilters": "-F 256 -F 512 -F 4 -F 2048",
        "ClipperSettings_nanopore": "",
        "ClipperSettings_illumina": "",
        "ClipperSettings_iontorrent": "",
        "Fastp_PhredScore_cutoff_illumina": 20,
        "Fastp_PhredScore_cutoff_nanopore": 7,
        "Fastp_PhredScore_cutoff_iontorrent": 20,
        "Fastp_WindowSize_illumina": 5,
        "Fastp_WindowSize_nanopore": 20,
        "Fastp_WindowSize_iontorrent": 15,
        "Fastp_MinReadLength": 100,
        "CleanAlign_AdditionalSettings_nanopore": "-E2,0 -O8,24 -A4 -B4",
        "CleanAlign_AdditionalSettings_illumina": "",
        "CleanAlign_AdditionalSettings_iontorrent": "",
    },
    "INFLUENZA": {
        "MatchRef_AdditionalAlignmentSettings": "--secondary=no",
        "MatchRef_AdditionalBamFilters": "-q 50",
        "RawAlign_AdditionalSettings": "--splice --frag=no",
        "BaseBAMFilters": "-F 256 -F 512 -F 4 -F 2048",
        "ClipperSettings_nanopore": "--exclude-spliced --spliced-length-threshold 50 --min-aligned-length 0.5",
        "ClipperSettings_illumina": "--exclude-spliced --spliced-length-threshold 50",
        "ClipperSettings_iontorrent": "--exclude-spliced --spliced-length-threshold 50 --min-aligned-length 0.5",
        "Fastp_PhredScore_cutoff_illumina": 20,
        "Fastp_PhredScore_cutoff_nanopore": 7,
        "Fastp_PhredScore_cutoff_iontorrent": 20,
        "Fastp_WindowSize_illumina": 5,
        "Fastp_WindowSize_nanopore": 20,
        "Fastp_WindowSize_iontorrent": 15,
        "Fastp_MinReadLength": 100,
        "CleanAlign_AdditionalSettings_nanopore": "",
        "CleanAlign_AdditionalSettings_illumina": "",
        "CleanAlign_AdditionalSettings_iontorrent": "",
    },
    "MEASLES": {
        "MatchRef_AdditionalAlignmentSettings": "--secondary=no",
        "MatchRef_AdditionalBamFilters": "-q 50",
        "RawAlign_AdditionalSettings": "",
        "BaseBAMFilters": "-F 256 -F 512 -F 4 -F 2048",
        "ClipperSettings_nanopore": "",
        "ClipperSettings_illumina": "",
        "ClipperSettings_iontorrent": "",
        "Fastp_PhredScore_cutoff_illumina": 20,
        "Fastp_PhredScore_cutoff_nanopore": 7,
        "Fastp_PhredScore_cutoff_iontorrent": 20,
        "Fastp_WindowSize_illumina": 5,
        "Fastp_WindowSize_nanopore": 20,
        "Fastp_WindowSize_iontorrent": 15,
        "Fastp_MinReadLength": 100,
        "CleanAlign_AdditionalSettings_nanopore": "",
        "CleanAlign_AdditionalSettings_illumina": "",
        "CleanAlign_AdditionalSettings_iontorrent": "",
    },
    "HPV": {
        "MatchRef_AdditionalAlignmentSettings": "--secondary=no",
        "MatchRef_AdditionalBamFilters": "-q 50",
        "RawAlign_AdditionalSettings": "",
        "BaseBAMFilters": "-F 256 -F 512 -F 4 -F 2048",
        "ClipperSettings_nanopore": "",
        "ClipperSettings_illumina": "",
        "ClipperSettings_iontorrent": "",
        "Fastp_PhredScore_cutoff_illumina": 20,
        "Fastp_PhredScore_cutoff_nanopore": 7,
        "Fastp_PhredScore_cutoff_iontorrent": 20,
        "Fastp_WindowSize_illumina": 5,
        "Fastp_WindowSize_nanopore": 20,
        "Fastp_WindowSize_iontorrent": 15,
        "Fastp_MinReadLength": 100,
        "CleanAlign_AdditionalSettings_nanopore": "",
        "CleanAlign_AdditionalSettings_illumina": "",
        "CleanAlign_AdditionalSettings_iontorrent": "",
    },
}


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


def get_preset_parameter(preset_name: str, parameter_name: str) -> str:
    """This function takes in a preset name and a parameter name, and returns the corresponding value for
    that parameter in the preset dictionary.

    Parameters
    ----------
    preset_name : str
        A string representing the name of a preset. Presets are pre-defined sets of values for certain
    parameters.
    parameter_name : str
        The name of the parameter that we want to retrieve from the preset.

    Returns
    -------
        A string value that corresponds to the parameter name of a given preset name. The value is
    retrieved from a dictionary called "presets".

    """
    return presets[preset_name][parameter_name]
