"""
Starting point of the ViroConstrictor pipeline and wrapper

Copyright © 2021 RIVM

https://github.com/RIVM-bioinformatics/ViroConstrictor
"""

# pylint: disable=C0103

import os
import sys

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

    ##> Check the default userprofile, make it if it doesn't exist
    conf = ReadConfig(os.path.expanduser("~/.ViroConstrictor_defaultprofile.ini"))

    flags, sampleinfo = ValidArgs(sys.argv[1:])

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
        workdir, inpath, start_path, conf, run_params, run_config, workflow_state,
    )

    if status is False:
        exit(1)
    exit(0)
