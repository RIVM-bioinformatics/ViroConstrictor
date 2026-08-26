import difflib
import json
import re
from typing import Any
from pathlib import Path
from enum import Enum

_HELPERS_DIR = Path(__file__).resolve().parent
aliases = json.loads((_HELPERS_DIR / "preset_aliases.json").read_text(encoding="utf-8"))
main_presets = json.loads((_HELPERS_DIR / "main_params.json").read_text(encoding="utf-8"))
matchref_presets = json.loads((_HELPERS_DIR / "matchref_params.json").read_text(encoding="utf-8"))


class Stage(Enum):
    STAGE_MAIN = "MAIN"
    STAGE_MR = "MATCHREF"

    @staticmethod
    def from_string(stage_str: str) -> "Stage":
        stage_str = stage_str.upper().strip()
        for stage in Stage:
            if stage.value == stage_str:
                return stage
        raise ValueError(f"Invalid stage identifier: {stage_str}")

    def __str__(self) -> str:
        return self.value


def _normalize_stage_identifier(stage_identifier: str | None) -> Stage:
    if stage_identifier:
        return Stage.from_string(stage_identifier)

    # backwards compatibility
    import inspect

    caller_globals = inspect.stack()[2][0].f_globals
    if "VC_STAGE" in caller_globals:
        return Stage.from_string(caller_globals["VC_STAGE"])
    raise ValueError("Stage identifier must be provided either as an argument or through the VC_STAGE global variable in the caller's scope.")


def match_preset_name(targetname: str, use_presets: bool) -> tuple[str, float]:
    """
    The function takes a target name and a boolean flag as input, and returns a tuple containing the
    best matching preset name and a score based on string similarity, or a default value if the flag is
    False or no match is found with a high enough similarity score.
    """
    if not use_presets:
        return "DEFAULT", 0.0

    query = re.sub(r"[^_a-zA-Z0-9/-]+", "", targetname).upper()
    if query == "DEFAULT":
        return "DEFAULT", 1.0

    aliases_list = [item for group in aliases.values() for item in group]
    if not aliases_list:
        return "DEFAULT", 0.0

    best = difflib.get_close_matches(query, aliases_list, cutoff=0.0, n=1)
    if not best:
        return "DEFAULT", 0.0

    best_match = best[0]
    score = difflib.SequenceMatcher(None, a=query, b=best_match).ratio()
    if score < 0.40:
        return "DEFAULT", 0.0

    for preset, values in aliases.items():
        if best_match in values:
            return preset, score
    return "DEFAULT", 0.0


def _search_presets(preset_name: str, task: str, platform: str, bin_or_flags: str) -> Any:
    preset = main_presets.get(preset_name, {}).get(task, {})
    result = preset.get(f"{bin_or_flags}_by_platform", {}).get(platform)
    if result is not None:
        return result
    # Try to get the result from the DEFAULT preset for the same platform
    default_preset = main_presets.get("DEFAULT", {}).get(task, {})
    result = default_preset.get(f"{bin_or_flags}_by_platform", {}).get(platform)
    if result is not None:
        return result
    raise KeyError(
        f"{bin_or_flags.capitalize()} not found. Preset: '{preset_name}', Platform: '{platform}', Task: '{task}', and no default {bin_or_flags} found for platform '{platform}'."
    )


def get_binary(stage_identifier: str, preset_name: str, task: str, platform: str) -> str:
    stage = _normalize_stage_identifier(stage_identifier)
    if stage == Stage.STAGE_MAIN:
        return _search_presets(preset_name, task, platform, "bin")
    elif stage == Stage.STAGE_MR:
        return _search_presets(preset_name, task, platform, "bin")
    else:
        raise ValueError(f"Invalid stage identifier: {stage_identifier}")


def get_flags(stage_identifier: str, preset_name: str, task: str, platform: str) -> str:
    stage = _normalize_stage_identifier(stage_identifier)
    if stage == Stage.STAGE_MAIN:
        return _search_presets(preset_name, task, platform, "flags")
    elif stage == Stage.STAGE_MR:
        return _search_presets(preset_name, task, platform, "flags")
    else:
        raise ValueError(f"Invalid stage identifier: {stage_identifier}")
