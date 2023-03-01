"""
Starting point of the ViroConstrictor pipeline and wrapper

Copyright © 2021 RIVM

https://github.com/RIVM-bioinformatics/ViroConstrictor
"""

# pylint: disable=C0103

import sys
from typing import Literal, NoReturn

import pandas as pd
import snakemake
from rich.progress import BarColumn, Progress, TextColumn

import ViroConstrictor.logging
from ViroConstrictor import __version__
from ViroConstrictor.logging import log
from ViroConstrictor.parser import CLIparser
from ViroConstrictor.runconfigs import GetSnakemakeRunDetails, WriteYaml
from ViroConstrictor.runreport import WriteReport
from ViroConstrictor.update import update


def get_preset_warning_list(
    sample_info_df: pd.DataFrame,
) -> tuple[list[str], list[str]]:
    """Takes a dataframe with sample information and returns a tuple of two lists of warnings

    Parameters
    ----------
    sample_info_df : pd.DataFrame
        pd.DataFrame

    Returns
    -------
        A list of warnings.

    """
    preset_fallback_warnings = preset_score_warnings = []
    for s in sample_info_df.itertuples():
        sample, preset, score, input_target = (
            s.SAMPLE,
            s.PRESET,
            s.PRESET_SCORE,
            s.VIRUS,
        )
        if score == 0.0:
            warn = f"""[red]Sample '[bold underline]{sample}[/bold underline]' was given the following information as an input-target: '[bold underline]{input_target}[/bold underline]'.
This information could however not be used to determine a preset. Because of this, the preset '[bold underline]{preset}[/bold underline]' was used instead.[/red]
[yellow]Please check your input-target for any significant misspellings, or consider using an alias or a different abbreviation for your input-target to check whether this resolves the issue.[/yellow]

It may be that your input-target does not yet have an associated preset in ViroConstrictor. 
If your suspect this to be the case, please open an issue on the ViroConstrictor GitHub page: [magenta underline]https://github.com/RIVM-bioinformatics/ViroConstrictor[/magenta underline]"""
            preset_fallback_warnings.append(warn)
            continue
        if score < 0.8:
            warn = f"""[red]Sample '[bold underline]{sample}[/bold underline]' was given the following information as an input-target: '[bold underline]{input_target}[/bold underline]'.
As a result, the preset '{preset}' was chosen. But this was done with less than 80% certainty.[/red]
[yellow]Certainty score: [bold]{score:.0%}[/bold][/yellow]
Please check the input-target and try again if a different preset is required."""
            preset_score_warnings.append(warn)
            continue
    return preset_fallback_warnings, preset_score_warnings


def show_preset_warnings(
    warnings: list[str], fallbacks: list[str], disabled: bool
) -> None:
    if warnings and not disabled:
        for w in warnings:
            log.warn(f"{w}")
    if fallbacks and not disabled:
        for w in fallbacks:
            log.warn(f"{w}")


def main() -> NoReturn:
    """
    ViroConstrictor starting point
    --> Fetch and parse arguments
    --> check validity
    --> Read (or write, if necessary) the user-config files
    --> Change working directories and make necessary local files for snakemake
    --> Run snakemake with appropriate settings
    """

    parsed_input = CLIparser(input_args=sys.argv[1:])

    preset_fallback_warnings, preset_score_warnings = get_preset_warning_list(
        parsed_input.samples_df
    )
    if not parsed_input.flags.skip_updates:
        update(sys.argv, parsed_input.user_config)

    snakemake_run_details = GetSnakemakeRunDetails(inputs_obj=parsed_input)

    log.info(f"{'='*25} [bold yellow] Starting Workflow [/bold yellow] {'='*25}")
    status: bool = False

    with Progress(
        TextColumn("[bold yellow]{task.description}[/bold yellow]"),
        BarColumn(bar_width=None),
        transient=True,
        expand=True,
    ) as progress:
        if parsed_input.user_config["COMPUTING"]["compmode"] == "local":
            progress.add_task(
                "Running ViroConstrictor workflow in local execution mode"
                if snakemake_run_details.snakemake_run_conf["dryrun"] is False
                else "Running Workflow in 'dryrun' mode",
                total=None,
            )
            status = snakemake.snakemake(
                snakefile=parsed_input.snakefile,
                workdir=parsed_input.workdir,
                cores=snakemake_run_details.snakemake_run_conf["cores"],
                use_conda=snakemake_run_details.snakemake_run_conf["use-conda"],
                conda_frontend="mamba",
                jobname=snakemake_run_details.snakemake_run_conf["jobname"],
                latency_wait=snakemake_run_details.snakemake_run_conf["latency-wait"],
                dryrun=snakemake_run_details.snakemake_run_conf["dryrun"],
                configfiles=[
                    WriteYaml(
                        snakemake_run_details.snakemake_run_parameters,
                        f"{parsed_input.workdir}/config/run_params.yaml",
                    )
                ],
                restart_times=3,
                keepgoing=True,
                quiet=["all"],  # type: ignore
                log_handler=[
                    ViroConstrictor.logging.snakemake_logger(
                        logfile=parsed_input.logfile
                    )
                ],
                printshellcmds=False,
            )
        if parsed_input.user_config["COMPUTING"]["compmode"] == "grid":
            progress.add_task(
                "Running ViroConstrictor workflow in grid execution mode"
                if snakemake_run_details.snakemake_run_conf["dryrun"] is False
                else "Running Workflow in 'dryrun' mode",
                total=None,
            )
            status = snakemake.snakemake(
                snakefile=parsed_input.snakefile,
                workdir=parsed_input.workdir,
                cores=snakemake_run_details.snakemake_run_conf["cores"],
                nodes=snakemake_run_details.snakemake_run_conf["cores"],
                use_conda=snakemake_run_details.snakemake_run_conf["use-conda"],
                conda_frontend="mamba",
                jobname=snakemake_run_details.snakemake_run_conf["jobname"],
                latency_wait=snakemake_run_details.snakemake_run_conf["latency-wait"],
                drmaa=snakemake_run_details.snakemake_run_conf["drmaa"],
                drmaa_log_dir=snakemake_run_details.snakemake_run_conf["drmaa-log-dir"],
                dryrun=snakemake_run_details.snakemake_run_conf["dryrun"],
                configfiles=[
                    WriteYaml(
                        snakemake_run_details.snakemake_run_parameters,
                        f"{parsed_input.workdir}/config/run_params.yaml",
                    )
                ],
                restart_times=3,
                keepgoing=True,
                quiet=["all"],  # type: ignore
                log_handler=[
                    ViroConstrictor.logging.snakemake_logger(
                        logfile=parsed_input.logfile
                    )
                ],
            )

    if snakemake_run_details.snakemake_run_conf["dryrun"] is False and status is True:
        snakemake.snakemake(
            snakefile=parsed_input.snakefile,
            workdir=parsed_input.workdir,
            report="results/snakemake_report.html",
            configfiles=[
                WriteYaml(
                    snakemake_run_details.snakemake_run_parameters,
                    f"{parsed_input.workdir}/config/run_params.yaml",
                )
            ],
            quiet=["all"],  # type: ignore
            log_handler=[
                ViroConstrictor.logging.snakemake_logger(logfile=parsed_input.logfile)
            ],
        )

    log.info(f"{'='*25} [bold yellow] Finished Workflow [/bold yellow] {'='*25}")
    workflow_state: Literal["Failed", "Success"] = (
        "Failed" if status is False else "Success"
    )

    WriteReport(
        parsed_input.workdir,
        parsed_input.input_path,
        parsed_input.exec_start_path,
        parsed_input.user_config,
        snakemake_run_details.snakemake_run_parameters,
        snakemake_run_details.snakemake_run_conf,
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
