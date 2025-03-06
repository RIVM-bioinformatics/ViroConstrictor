import datetime
import logging
import os
import pathlib
import re
from typing import Any, Callable

from rich.color import ANSI_COLOR_NAMES
from rich.default_styles import DEFAULT_STYLES
from rich.highlighter import NullHighlighter
from rich.logging import RichHandler

richstyles = list(set(list(DEFAULT_STYLES.keys()) + list(ANSI_COLOR_NAMES.keys())))


class StripBracketsFilter(logging.Filter):
    """
    A class used to strip rich-style markup from log messages before they are written to a log file.
    """

    def filter(self, record):
        pattern = r"\[/?\b(%s)\b[^\]]*\]" % "|".join(richstyles)
        record.msg = re.sub(pattern, "", record.msg)
        record.msg = record.msg.replace("\n", "\n\t\t\t\t\t\t\t")
        return True


log = logging.getLogger("ViroConstrictor")
log.setLevel("INFO")


def setup_logger(workdir: str) -> str:
    """This function sets up a logger for ViroConstrictor and returns the path to the log file.

    Parameters
    ----------
    workdir : str
        The `workdir` parameter is a string representing the directory path where the log file will be
    saved.

    Returns
    -------
        A string which is the absolute path of the log file created by the function.

    """
    logging.getLogger("asyncio").setLevel(logging.CRITICAL)
    logging.getLogger("snakemake").setLevel(logging.CRITICAL)
    logging.getLogger("smart_open").setLevel(logging.CRITICAL)
    logging.getLogger("urllib3").setLevel(logging.CRITICAL)
    logging.getLogger("fpdf").setLevel(logging.CRITICAL)

    if not os.path.exists(workdir):
        os.makedirs(workdir)
    logfile = f"{os.path.abspath(workdir)}/ViroConstrictor_{datetime.datetime.now().isoformat()}.log"

    shell_handler = RichHandler(
        show_path=False,
        omit_repeated_times=False,
        markup=True,
        highlighter=NullHighlighter(),
        rich_tracebacks=True,
    )
    file_handler = logging.FileHandler(logfile)

    shell_handler.setLevel("INFO")
    file_handler.setLevel("INFO")

    format_shell = "%(message)s"
    format_file = "%(asctime)s\t%(levelname)s\t%(message)s"

    shell_formatter = logging.Formatter(format_shell, "[%d/%m/%y %H:%M:%S]")
    file_formatter = logging.Formatter(format_file, "[%d/%m/%y %H:%M:%S]")

    shell_handler.setFormatter(shell_formatter)
    file_handler.addFilter(StripBracketsFilter())
    file_handler.setFormatter(file_formatter)

    log.addHandler(shell_handler)
    log.addHandler(file_handler)

    return logfile


def BaseLogMessage(msg: dict[str, Any]) -> None:
    if logmessage := msg.get("msg"):
        log.info(f"{logmessage}")


def CondaEnvInstallerPreamble(msg: dict[str, Any]) -> None:
    if logmessage := msg.get("msg"):
        log.info(f"Creating conda environment: [yellow]{logmessage}[/yellow]")


def CondaEnvInstallSuccess(msg: dict[str, Any]) -> None:
    log.info("[green]Conda environment created![/green]")


def BaseSnakemakeAbortMessage(msg: dict[str, Any]) -> None:
    if logmessage := msg.get("msg"):
        log.error(f"[red]{logmessage}[/red]")


def ColorizeLogMessagePath(msg: dict[str, Any]) -> None:
    """This function colorizes any file paths found in a log message with magenta color.

    Parameters
    ----------
    msg : dict[str, Any]
        The `msg` parameter is a dictionary that contains information about a log message.

    """
    if logmessage := msg.get("msg"):
        msgparts = logmessage.split(" ")
        for word in msgparts:
            word = word.replace("...", "") if "..." in word else word.rstrip(".")
            if pathlib.Path(word).exists():
                logmessage = logmessage.replace(word, f"[magenta]{word}[/magenta]")
        log.info(f"{logmessage}")


def LogFileOverride(msg: dict[str, Any], logfile) -> None:
    """This function logs a message with the name of a logfile.

    Parameters
    ----------
    msg : dict[str, Any]
        The `msg` parameter is a dictionary that contains the original logging dictionary from snakemake. It is expected to have
    a key "msg" which contains the original log message retrieved from snakemake. This serves as a trigger for the function to run.
    logfile
        The `logfile` parameter is a variable that represents the file path or name of the log file that
    will be used to store the log messages.

    """
    if msg.get("msg") is not None:
        log.info(f"Complete log: [magenta]{logfile}[/magenta]")


def SubmitDRMAAMessage(msg: dict[str, Any]) -> None:
    """The function extracts job identifiers from a message dictionary and logs a message indicating that a
    job has been submitted to a DRMAA cluster using an external job ID.

    Parameters
    ----------
    msg : dict[str, Any]
        The parameter `msg` is a dictionary that contains a log message.

    """
    if logmessage := msg.get("msg"):
        logmessage_items = logmessage.rstrip(".").split(" ")
        identifiers = []
        for item in logmessage_items:
            try:
                identifiers.append(int(item))
            except Exception:
                continue
        log.info(
            f"Submitted JobID [cyan]{identifiers[0]}[/cyan] to the DRMAA cluster using external JobID [cyan]{identifiers[1]}[/cyan]."
        )


def HandleJobInfoMessage(msg: dict[str, Any]) -> None:
    """The function logs information about a job being executed, including the process name, sample, input
    target, reference ID, job ID, and log file.

    Parameters
    ----------
    msg : dict[str, Any]
        The `msg` parameter is a dictionary that contains information about a job to be executed.

    """
    if not (processname := msg.get("name")):
        return
    wildcards = msg.get("wildcards")
    sample = wildcards.get("sample") if wildcards else None
    input_target = wildcards.get("Virus") if wildcards else None
    refid = wildcards.get("RefID") if wildcards else None
    logfile = msg.get("log")
    logfile = pathlib.Path.absolute(pathlib.Path(logfile[0])) if logfile else None
    if msg.get("local"):
        log.info(
            f"Executing localjob [green underline]{processname}[/green underline] for sample [blue]{sample}[/blue] with target [blue]{input_target}[/blue] and reference-id [blue]{refid}[/blue]\nJob is using jobID [cyan]{msg.get('jobid')}[/cyan], logging output will be written to [magenta]{logfile}[/magenta]"
        )
    else:
        log.info(
            f"Executing job [green underline]{processname}[/green underline] for sample [blue]{sample}[/blue] with target [blue]{input_target}[/blue] and reference-id [blue]{refid}[/blue]\nJob is using jobID [cyan]{msg.get('jobid')}[/cyan], logging output will be written to [magenta]{logfile}[/magenta]"
        )


def HandleJobErrorMessage(msg: dict[str, Any]) -> None:
    """This function handles error messages for a job and logs relevant information.

    Parameters
    ----------
    msg : dict[str, Any]
        The `msg` parameter is a dictionary containing information about a job error message.

    """
    if not (processname := msg.get("name")):
        return
    wildcards = msg.get("wildcards")
    sample = wildcards.get("sample") if wildcards else None
    input_target = wildcards.get("Virus") if wildcards else None
    refid = wildcards.get("RefID") if wildcards else None
    logfile = msg.get("log")
    logfile = pathlib.Path.absolute(pathlib.Path(logfile[0])) if logfile else None
    shellcmd = str(msg.get("shellcmd")).strip()
    if outputfiles := msg.get("output"):
        outputfiles_list = " ".join(list(outputfiles))
    else:
        outputfiles_list = "None"
    log.error(
        f"Job [red underline]{processname}[/red underline] for sample [blue]{sample}[/blue] with target [blue]{input_target}[/blue] and reference-id [blue]{refid}[/blue] failed!\nJob is using jobID [cyan]{msg.get('jobid')}[/cyan], logging output will be written to [magenta]{logfile}[/magenta]\nThe expected output file(s) are as follows:\n[red]{outputfiles_list}[/red]\nThe following shell command was issued:\n[red]{shellcmd}[/red]"
    )


def print_jobstatistics_logmessage(msg: dict) -> None:
    if logmessage := msg.get("msg"):
        logmessage = logmessage.split("\n", 1)[1]
        log.info(f"Job statistics:\n[yellow]{logmessage}[/yellow]")


logmessage_strings_info: dict[str, Callable] = {
    "Activating conda environment": ColorizeLogMessagePath,
    "Activating singularity image": ColorizeLogMessagePath,
    "Building DAG of jobs": BaseLogMessage,
    "Creating conda environment": ColorizeLogMessagePath,
    "Removing incomplete Conda environment": ColorizeLogMessagePath,
    "Downloading and installing remote packages": CondaEnvInstallerPreamble,
    "Environment for": CondaEnvInstallSuccess,
    "Cancelling snakemake on user request": BaseSnakemakeAbortMessage,
    "Select jobs to execute...": BaseLogMessage,
    "Creating report...": BaseLogMessage,
    "Downloading resources and rendering HTML.": BaseLogMessage,
    "Report created: results/snakemake_report.html.": ColorizeLogMessagePath,
    "Terminating processes on user request, this might take some time.": BaseSnakemakeAbortMessage,
    "Removing temporary output": ColorizeLogMessagePath,
    "will be created.": ColorizeLogMessagePath,
    "Submitted DRMAA job": SubmitDRMAAMessage,
    "Nothing to be done (all requested files are present and up to date).": BaseLogMessage,
}
logmessage_suppressed_strings_warning: list[str] = [
    "Your conda installation is not configured to use strict channel priorities.",
    "Failed to download resource needed for report:"
]


def snakemake_logger(logfile: str) -> object:
    """`snakemake_logger` is a wrapper function that returns a log handler that processes different types of log messages
    and outputs them accordingly.
    This function has the purpose of allowing the `logfile` variable to be passed to the underlying log_handler.

    Parameters
    ----------
    logfile : str
        The path and filename of the log file where the log messages will be written.

    Returns
    -------
        The function `snakemake_logger` returns a function object `log_handler`.

    """

    def log_handler(msg: dict) -> None:
        loglevel = msg.get("level")
        logmessage = msg.get("msg")

        if loglevel in ["dag_debug", "debug", "shellcmd"]:
            return None

        if logmessage is not None and any(
            x in logmessage for x in logmessage_suppressed_strings_warning
        ):
            return None

        if logmessage is not None and any(
            x in logmessage for x in list(logmessage_strings_info.keys())
        ):
            for key in logmessage_strings_info.keys():
                if key in logmessage:
                    logmessage_strings_info[key](msg)

        elif logmessage is not None and "Complete log:" in logmessage:
            LogFileOverride(msg, logfile)
        elif loglevel == "run_info":
            print_jobstatistics_logmessage(msg)
        elif loglevel == "job_info":
            HandleJobInfoMessage(msg)
        elif loglevel == "job_error":
            HandleJobErrorMessage(msg)
        elif loglevel == "warning":
            log.warning(f"{logmessage}")
        elif loglevel == "error":
            log.error(f"{logmessage}")

    return log_handler
