"""
Starting point of the ViroConstrictor pipeline and wrapper

Copyright © 2021 RIVM

https://github.com/RIVM-bioinformatics/ViroConstrictor
"""

# pylint: disable=C0103

import os
import sys

import pandas as pd
import snakemake
import yaml

from ViroConstrictor import __version__
from ViroConstrictor.functions import color
from ViroConstrictor.parser import ValidArgs
from ViroConstrictor.runconfigs import SnakemakeConfig, SnakemakeParams, WriteYaml
from ViroConstrictor.runreport import WriteReport
from ViroConstrictor.update import update
from ViroConstrictor.userprofile import ReadConfig
from ViroConstrictor.validatefasta import CheckReferenceFile

yaml.warnings({"YAMLLoadWarning": False})


def CheckSampleProperties(sampleinfo):
    """This function checks that the reference fasta file exists and that it is a valid fasta file

    Parameters
    ----------
    sampleinfo
        A dictionary of dictionaries. The outer dictionary is keyed by sample name, and the inner
        dictionary should contain the "REFERENCE" key containing a path to the reference fasta file.

    """
    for sample in sampleinfo:
        if not os.path.isfile(sampleinfo.get(sample).get("REFERENCE")):
            raise FileNotFoundError(
                f"\n{color.RED + color.BOLD}The given reference fasta file for sample '{sample}' does not exist. Please check the reference fasta and try again. Exiting...{color.END}\n"
            )

    reference_files = {sampleinfo.get(sample).get("REFERENCE") for sample in sampleinfo}
    for f in reference_files:
        CheckReferenceFile(f)


def main():
    """
    ViroConstrictor starting point
    --> Fetch and parse arguments
    --> check validity
    --> Read (or write, if necessary) the user-config files
    --> Change working directories and make necessary local files for snakemake
    --> Run snakemake with appropriate settings
    """

    flags, sampleinfo, samples_df = ValidArgs(sys.argv[1:])
    samples_df = samples_df.reset_index(drop=False).rename(columns={"index": "SAMPLE"})
    sampleinfo_df = (
        pd.DataFrame.from_dict(sampleinfo, orient="index")
        .reset_index(drop=False)
        .rename(columns={"index": "SAMPLE"})
    )

    ##> Check the default userprofile, make it if it doesn't exist
    conf = ReadConfig(os.path.expanduser("~/.ViroConstrictor_defaultprofile.ini"))

    preset_fallback_warnings = []
    preset_score_warnings = []
    for s in sampleinfo_df.itertuples():
        sample, preset, score, input_target = (
            s.SAMPLE,
            s.PRESET,
            s.PRESET_SCORE,
            s.VIRUS,
        )
        if score == 0.0:
            warn = Warning(
                f"""Sample '{sample}' was given the following information as an input-target: '{input_target}'.
This information could however not be used to determine a preset. Because of this, the preset '{preset}' was used instead.
Please check your input-target for any significant misspellings, or consider using an alias or a different abbreviation for your input-target to check whether this resolves the issue.

It may be that your input-target does not yet have an associated preset in ViroConstrictor. 
If your suspect this to be the case, please open an issue on the ViroConstrictor GitHub page: https://github.com/RIVM-bioinformatics/ViroConstrictor"""
            )
            preset_fallback_warnings.append(warn)
            continue
        if score < 0.8:
            warn = Warning(
                f"""Sample '{sample}' was given the following information as an input-target: '{input_target}'.
As a result, the preset '{preset}' was chosen. But this was done with less than 80% certainty.
Certainty score: {score:.0%}
Please check the input-target and try again if a different preset is required."""
            )
            preset_score_warnings.append(warn)
            continue

    if not flags.skip_updates:
        update(sys.argv, conf)

    inpath = os.path.abspath(flags.input)
    start_path = os.getcwd()

    outpath = os.path.abspath(flags.output)

    exec_folder = os.path.abspath(os.path.dirname(__file__))

    Snakefile = os.path.join(exec_folder, "workflow", "workflow.smk")

    CheckSampleProperties(sampleinfo)  # raises errors if stuff is not right

    ##@ check if the output dir exists, create if not
    ##@ change the working directory
    if not os.path.exists(outpath):
        os.makedirs(outpath)

    if os.getcwd() != outpath:
        os.chdir(outpath)
    workdir = outpath

    samplesheet = WriteYaml(sampleinfo, f"{workdir}/samplesheet.yaml")
    run_config = SnakemakeConfig(conf, flags.threads, flags.dryrun)
    run_params = SnakemakeParams(
        conf,
        flags.threads,
        sampleinfo,
        flags.platform,
        samplesheet,
        flags.amplicon_type,
    )

    if conf["COMPUTING"]["compmode"] == "local":
        status = snakemake.snakemake(
            Snakefile,
            workdir=workdir,
            cores=run_config["cores"],
            use_conda=run_config["use-conda"],
            conda_frontend="mamba",
            jobname=run_config["jobname"],
            latency_wait=run_config["latency-wait"],
            dryrun=run_config["dryrun"],
            configfiles=[WriteYaml(run_params, f"{workdir}/config/run_params.yaml")],
            restart_times=3,
            keepgoing=True,
        )
    if conf["COMPUTING"]["compmode"] == "grid":
        status = snakemake.snakemake(
            Snakefile,
            workdir=workdir,
            cores=run_config["cores"],
            nodes=run_config["cores"],
            use_conda=run_config["use-conda"],
            conda_frontend="mamba",
            jobname=run_config["jobname"],
            latency_wait=run_config["latency-wait"],
            drmaa=run_config["drmaa"],
            drmaa_log_dir=run_config["drmaa-log-dir"],
            dryrun=run_config["dryrun"],
            configfiles=[WriteYaml(run_params, f"{workdir}/config/run_params.yaml")],
            restart_times=3,
            keepgoing=True,
        )

    if run_config["dryrun"] is False and status is True:
        snakemake.snakemake(
            Snakefile,
            workdir=workdir,
            report="results/snakemake_report.html",
            configfiles=[WriteYaml(run_params, f"{workdir}/config/run_params.yaml")],
            quiet=True,
        )

    workflow_state = "Failed" if status is False else "Success"

    WriteReport(
        workdir,
        inpath,
        start_path,
        conf,
        run_params,
        run_config,
        workflow_state,
    )

    if preset_score_warnings and not flags.disable_presets:
        for w in preset_score_warnings:
            print(f"{color.YELLOW}Warning: {w}{color.END}\n")
    if preset_fallback_warnings and not flags.disable_presets:
        for w in preset_fallback_warnings:
            print(f"{color.YELLOW + color.BOLD}Warning: {w}{color.END}\n")

    if status is False:
        exit(1)
    exit(0)
