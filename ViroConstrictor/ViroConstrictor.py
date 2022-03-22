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
from ViroConstrictor.runconfigs import LoadConf, WriteConfigs
from ViroConstrictor.runreport import WriteReport
from ViroConstrictor.samplesheet import WriteSampleSheet
from ViroConstrictor.update import update
from ViroConstrictor.userprofile import ReadConfig
from ViroConstrictor.validatefasta import IsValidFasta

yaml.warnings({"YAMLLoadWarning": False})


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
    refpath = os.path.abspath(flags.reference)

    if flags.primers != "NONE":
        primpath = os.path.abspath(flags.primers)
    else:
        primpath = "NONE"

    if flags.features != "NONE":
        featspath = os.path.abspath(flags.features)
    else:
        featspath = "NONE"

    outpath = os.path.abspath(flags.output)

    exec_folder = os.path.abspath(os.path.dirname(__file__))

    Snakefile = os.path.join(exec_folder, "workflow", "workflow.smk")

    ##@ check if the input directory contains valid files
    if CheckInputFiles(inpath) is False:
        print(
            f"""
{color.RED + color.BOLD}"{inpath}" does not contain any valid FastQ files.{color.END}
Please check the input directory. Exiting...
            """
        )
        sys.exit(-1)
    else:
        print(
            f"""
{color.GREEN}Valid input files were found in the input directory{color.END} ('{inpath}')
            """
        )

    if IsValidFasta(primpath) is False:
        print(
            f"""
{color.RED + color.BOLD}
The given fasta with primer sequences contains illegal characters in its sequences.
{color.END}
Please check the primer fasta and try again. Exiting...
            """
        )
        sys.exit(1)

    if IsValidFasta(refpath) is False:
        print(
            f"""
{color.RED + color.BOLD}
The given fasta with the reference sequence contains illegal characters in its sequences.
{color.END}
Please check the reference fasta and try again. Exiting...
            """
        )
        sys.exit(1)

    ##@ check if the output dir exists, create if not
    ##@ change the working directory
    if not os.path.exists(outpath):
        os.makedirs(outpath)

    # copy_tree(os.path.join(here, 'envs'), os.path.join(outpath, 'envs'))

    if not os.getcwd() == outpath:
        os.chdir(outpath)
    workdir = outpath

    samplesheet = WriteSampleSheet(inpath, flags.platform)
    snakeparams, snakeconfig = WriteConfigs(
        conf,
        flags.threads,
        os.getcwd(),
        flags.platform,
        refpath,
        primpath,
        featspath,
        samplesheet,
        flags.amplicon_type,
        flags.primer_mismatch_rate,
        flags.dryrun,
    )

    parsedconfig = LoadConf(snakeconfig)

    if conf["COMPUTING"]["compmode"] == "local":
        status = snakemake.snakemake(
            Snakefile,
            workdir=workdir,
            cores=parsedconfig["cores"],
            use_conda=parsedconfig["use-conda"],
            conda_frontend="mamba",
            jobname=parsedconfig["jobname"],
            latency_wait=parsedconfig["latency-wait"],
            dryrun=parsedconfig["dryrun"],
            configfiles=[snakeparams],
            restart_times=3,
        )
    if conf["COMPUTING"]["compmode"] == "grid":
        status = snakemake.snakemake(
            Snakefile,
            workdir=workdir,
            cores=parsedconfig["cores"],
            nodes=parsedconfig["cores"],
            use_conda=parsedconfig["use-conda"],
            conda_frontend="mamba",
            jobname=parsedconfig["jobname"],
            latency_wait=parsedconfig["latency-wait"],
            drmaa=parsedconfig["drmaa"],
            drmaa_log_dir=parsedconfig["drmaa-log-dir"],
            dryrun=parsedconfig["dryrun"],
            configfiles=[snakeparams],
            restart_times=3,
        )

    if parsedconfig["dryrun"] is False and status is True:
        snakemake.snakemake(
            Snakefile,
            workdir=workdir,
            report="results/snakemake_report.html",
            configfiles=[snakeparams],
            quiet=True,
        )

    if status is False:
        workflow_state = "Failed"
    else:
        workflow_state = "Success"

    WriteReport(
        workdir,
        inpath,
        start_path,
        conf,
        LoadConf(snakeparams),
        LoadConf(snakeconfig),
        workflow_state,
    )

    if status is True:
        exit(0)
    else:
        exit(1)
