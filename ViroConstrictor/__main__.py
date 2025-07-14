"""
Starting point of the ViroConstrictor pipeline and wrapper

Copyright © 2021 RIVM

https://github.com/RIVM-bioinformatics/ViroConstrictor
"""

# pylint: disable=C0103

import logging
import sys
from itertools import zip_longest
from typing import Literal, NoReturn

import pandas as pd

from ViroConstrictor import __version__
from ViroConstrictor.logging import log
from ViroConstrictor.match_ref import process_match_ref
from ViroConstrictor.parser import CLIparser
from ViroConstrictor.runreport import WriteReport
from ViroConstrictor.update import update
from ViroConstrictor.workflow_executor import run_snakemake_workflow


def get_preset_warning_list(
    sample_info_df: pd.DataFrame,
) -> tuple[list[str], list[str]]:
    """This function generates warning messages for cases where the preset used in ViroConstrictor was
    determined with low certainty or could not be determined at all.

    Parameters
    ----------
    sample_info_df : pd.DataFrame
        `sample_info_df` is a pandas DataFrame containing information about the samples being analyzed,
    including the virus being targeted, the preset used for analysis, and the preset score (a measure of
    how well the preset matches the input-target).

    Returns
    -------
        a tuple containing two lists of strings: `preset_fallback_warnings` and `preset_score_warnings`.

    """
    preset_fallback_warnings = []
    preset_score_warnings = []

    p_scorewarning_df = sample_info_df.loc[
        (sample_info_df["PRESET_SCORE"] < 0.8) & (sample_info_df["PRESET_SCORE"] > 0.0)
    ]
    for _input, _preset in zip(
        list(set(p_scorewarning_df["VIRUS"].tolist())),
        list(set(p_scorewarning_df["PRESET"].tolist())),
    ):
        filtered_df = p_scorewarning_df.loc[
            (p_scorewarning_df["VIRUS"] == _input)
            & (p_scorewarning_df["PRESET"] == _preset)
        ]
        samples = [" * " + x + "\n" for x in filtered_df["SAMPLE"].tolist()]
        score = filtered_df["PRESET_SCORE"].tolist()[0]

        warn = f"""[red]The following information was given as an input-target: '[bold underline]{_input}[/bold underline]'.
As a result, the preset '[bold underline]{_preset}[/bold underline]' was chosen. But this was done with less than 80% certainty.[/red]
[yellow]Certainty score: [bold]{score:.0%}[/bold][/yellow]
Please check the input-target and try again if a different preset is required, or use the [italic]"--disable-presets"[/italic] flag in order to always use the default preset.
This applies to the following samples:\n{''.join(samples)}"""
        preset_score_warnings.append(warn)

    # check if the preset score is larger or equal than 0.0 and smaller than 0.000001 (1e-6)
    # We do this because the preset score is a float and we want to check if it is within a certain range as floating point equality checks are not reliable
    p_fallbackwarning_df = sample_info_df.loc[
        (sample_info_df["PRESET_SCORE"] >= 0.0)
        & (sample_info_df["PRESET_SCORE"] < 1e-6)
    ]

    targets, presets = (
        (
            list(x)
            for x in (
                zip(
                    *zip_longest(
                        list(set(p_fallbackwarning_df["VIRUS"].tolist())),
                        list(set(p_fallbackwarning_df["PRESET"].tolist())),
                        fillvalue="DEFAULT",
                    )
                )
            )
        )
        if p_fallbackwarning_df.shape[0] > 0
        else ([], [])
    )
    for _input, _preset in zip(targets, presets):
        filtered_df = p_fallbackwarning_df.loc[
            (p_fallbackwarning_df["VIRUS"] == _input)
            & (p_fallbackwarning_df["PRESET"] == _preset)
        ]
        samples = [" * " + x + "\n" for x in filtered_df["SAMPLE"].tolist()]

        warn = f"""[red]The following information was given as an input-target: '[bold underline]{_input}[/bold underline]'.
This information could however not be used to determine a preset. Because of this, the preset '[bold underline]{_preset}[/bold underline]' was used instead.[/red]
[yellow]Please check your input-target for any significant misspellings, or consider using an alias or a different abbreviation for your input-target to check whether this resolves the issue.[/yellow]
This applies to the following samples:\n{''.join(samples)}
It may also be possible that your input-target does not yet have an associated preset in ViroConstrictor. Please see the documentation for more information regarding the currently available presets."""
        preset_fallback_warnings.append(warn)
    return preset_fallback_warnings, preset_score_warnings


def show_preset_warnings(
    warnings: list[str], fallbacks: list[str], disabled: bool
) -> None:
    """This function logs warning and fallback messages if they exist and if the disabled flag is not set.

    Parameters
    ----------
    warnings : list[str]
        A list of warning messages to be displayed.
    fallbacks : list[str]
        The `fallbacks` parameter is a list of strings that contains warnings about preset-fallback behavior.
    disabled : bool
        The "disabled" parameter is a boolean flag that indicates whether or not warning and fallback messages
    should be displayed. If it is set to True, then no warnings or fallbacks will be shown. If it is set
    to False, then warnings and fallbacks will be shown if there are any.

    """
    if warnings and not disabled:
        for w in warnings:
            logging.warning(f"{w}")
    if fallbacks and not disabled:
        for w in fallbacks:
            logging.warning(f"{w}")


def main(args: list[str] | None = None, settings: str | None = None) -> NoReturn:
    """
    ViroConstrictor starting point
    --> Fetch and parse arguments
    --> check validity
    --> Read (or write, if necessary) the user-config files
    --> Change working directories and make necessary local files for snakemake
    --> Run snakemake with appropriate settings
    """
    if args is None:
        args = sys.argv[1:]

    if settings is None:
        settings = "~/.ViroConstrictor_defaultprofile.ini"
    parsed_input = CLIparser(input_args=args, settings_path=settings)

    preset_fallback_warnings, preset_score_warnings = get_preset_warning_list(
        parsed_input.samples_df
    )
    if not parsed_input.flags.skip_updates:
        update(sys.argv, parsed_input.user_config)

    # check if there's a value in the column 'MATCH-REF' set to True in the parsed_input.samples_df dataframe, if so, process the match-ref, else skip
    if parsed_input.samples_df["MATCH-REF"].any():
        parsed_input = process_match_ref(parsed_input, scheduler=parsed_input.scheduler)

    log.info(f"{'='*20} [bold yellow] Starting Main Workflow [/bold yellow] {'='*20}")

    status: bool = False

    status, used_workflow_config = run_snakemake_workflow(
        inputs_obj=parsed_input, stage="MAIN", scheduler=parsed_input.scheduler
    )

    # if used_workflow_config.output_settings.dryrun is False and status is True:
    #    #     TODO: add back the functionality to generate a snakemake report upon completion of the workflow

    workflow_state: Literal["Failed", "Success"] = (
        "Failed" if status is False else "Success"
    )

    WriteReport(
        parsed_input.workdir,
        parsed_input.input_path,
        parsed_input.exec_start_path,
        parsed_input.user_config,
        used_workflow_config.resource_settings,
        used_workflow_config.output_settings,
        parsed_input,
        workflow_state,
    )

    show_preset_warnings(
        preset_score_warnings,
        preset_fallback_warnings,
        parsed_input.flags.disable_presets,
    )

    if status is False:
        exit(1)
    exit(0)
