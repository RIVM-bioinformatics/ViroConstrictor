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

from ViroConstrictor import __prog__

richstyles = list(set(list(DEFAULT_STYLES.keys()) + list(ANSI_COLOR_NAMES.keys())))


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
    global logfile
    logfile = f"{os.path.abspath(workdir)}/ViroConstrictor_{datetime.datetime.now().isoformat()}.log"

    # Initialize the ViroConstrictorBaseLogHandler to set up the log handlers in a way that both the python application and snakemake emit to the same handlers.
    ViroConstrictorBaseLogHandler(logfile_path=logfile)

    return logfile


class StripBracketsFilter(logging.Filter):
    """
    A class used to strip rich-style markup from log messages before they are written to a log file.
    """

    def filter(self, record):
        pattern = r"\[/?\b(%s)\b[^\]]*\]" % "|".join(richstyles)
        record.msg = re.sub(pattern, "", record.msg)
        record.msg = record.msg.replace("\n", "\n\t\t\t\t\t\t\t")
        return record


class ViroConstrictorBaseLogHandler(logging.Handler):
    """
    A custom log handler for ViroConstrictor that extends the standard logging.Handler.
    This handler provides enhanced logging capabilities, including:
    - Simultaneous console and file output.
    - Rich formatting for console output using the `rich` library.
    - Filtering of log messages to remove brackets from file output.
    - Customizable formatting for file output.
    - Suppression of specific log events and messages.
    - Special handling for job-related log events (start, info, error).

    Attributes
    ----------
    _strip_brackets_filter_instance : StripBracketsFilter
        A static instance of the StripBracketsFilter used to remove brackets from log messages written to file.
        This is instantiated once and shared across all instances of this handler.
    console_handler : RichHandler
        An instance of RichHandler used for console output, configured with specific settings
        such as disabling path display, omitting repeated times, enabling markup, and using rich tracebacks.
    shell_handler : RichHandler
        An instance of RichHandler used for shell output, configured to avoid colorization in the console.
    instance_file_handler : logging.FileHandler
        A FileHandler instance specific to this log handler, responsible for writing log messages to a file.
        It includes a filter to strip brackets and uses a custom formatter.

    Methods
    --------
    emit(record)
        Emits a log record to both the console and a file, applying rich formatting to the console output
        and stripping brackets from the file output.
    close()
        Closes both the console and file handlers, ensuring that all resources are released.
    _dispatch_log_record_formatter(record)
    """

    # Instantiate the filter once to be shared by all instances or used per instance
    _strip_brackets_filter_instance: StripBracketsFilter = StripBracketsFilter()

    def __init__(self, logfile_path: str, *args, **kwargs):
        # Initialize the base Handler
        super().__init__()

        # This handler will manage its own console handler internally
        self.console_handler = RichHandler(
            show_path=False,
            omit_repeated_times=False,
            markup=True,
            rich_tracebacks=True,
            highlighter=NullHighlighter(),
            *args,
            **kwargs,
        )

        # clear all existing log handlers so we have a clean slate every time we initialize this handler.
        for handler in log.handlers[:]:
            log.removeHandler(handler)

        # Configure the console handler to use rich formatting
        self.shell_handler = RichHandler(
            show_path=False,
            omit_repeated_times=False,
            markup=True,
            rich_tracebacks=True,
            highlighter=NullHighlighter(),  # Use NullHighlighter to avoid colorization in the console
        )
        log.addHandler(self.shell_handler)

        # Setup a FileHandler specific to this instance for writing to the log file
        self.instance_file_handler = logging.FileHandler(logfile_path)
        self.instance_file_handler.addFilter(
            self._strip_brackets_filter_instance
        )  # Add the filter to the file handler

        # Define formatter for the file log
        format_file = "%(asctime)s\t%(levelname)s\t%(message)s"
        file_formatter = logging.Formatter(
            fmt=format_file, datefmt="[%d/%m/%y %H:%M:%S]"
        )
        self.instance_file_handler.setFormatter(file_formatter)

        log.addHandler(self.instance_file_handler)

    def emit(self, record: logging.LogRecord) -> None:
        """
        Emit a log record.
        This method is responsible for handling both console and file output of log records.
        It first processes the record for console output using the RichHandler's capabilities.
        Then, it creates a copy of the record, applies a filter to strip brackets from the message,
        and emits the modified record to the file.
        Parameters
        ----------
        record : logging.LogRecord
            The log record to be emitted.
        Returns
        -------
        None
        Raises
        ------
        Exception
            If any error occurs during console or file output, the error is handled by calling `self.handleError(record)`.
        Notes
        -----
        - The console output is handled by the internal `console_handler`.
        - The file output involves creating a copy of the record to avoid modifications affecting other handlers.
        - The `StripBracketsFilter` is applied to the file record to modify the message in place.
        - The file output is handled by the `instance_file_handler`.
        """

        # take the provided log record and format it for console output with rich formatting.
        processed_record = self._dispatch_log_record_formatter(record)
        if processed_record is None:
            return

        try:
            # Use the internal console_handler to emit to the console
            self.console_handler.emit(processed_record)
        except Exception:
            self.handleError(record)

        # 2. Handle file output
        # If self.emit() is called, the record should also be written to the file.
        try:
            # Create a copy of the record for file processing.
            # This is to prevent modifications from affecting other handlers or console output.
            file_record = logging.makeLogRecord(record.__dict__)

            # Get the fully resolved message (msg % args).
            # This is important so the filter operates on the final string.
            message_for_file = file_record.getMessage()

            # Place the fully resolved message into file_record.msg
            # and clear args, so the filter operates on the complete string.
            file_record.msg = message_for_file
            file_record.args = ()

            # Apply the StripBracketsFilter. This modifies file_record.msg in place.
            # This is necessary to ensure that the message written to the file does not contain rich-style markup.
            self._strip_brackets_filter_instance.filter(file_record)

            # Emit the modified record using the file handler.
            # FileHandler.emit() itself does not re-check the level.
            self.instance_file_handler.emit(file_record)
        except Exception:
            self.handleError(record)

    def close(self) -> None:
        try:
            # Close the internal console handler
            self.console_handler.close()
        finally:
            # Ensure the file handler is also closed
            if self.instance_file_handler:
                self.instance_file_handler.close()
        super().close()

    def _dispatch_log_record_formatter(
        self, record: logging.LogRecord
    ) -> logging.LogRecord | None:
        """
        Formats a log record based on its content, applying specific transformations or filters.
        This method inspects the log record's level, message, event, and associated log file,
        and applies transformations such as colorizing specific messages, suppressing certain events or messages,
        or handling specific events like job start, info, or error.

        Parameters
        ----------
        record : logging.LogRecord
            The log record to be formatted.

        Returns
        -------
        logging.LogRecord | None
            The formatted log record, or None if the record should be suppressed.

        Raises
        ------
        KeyError
            If a required key is missing from the log record.
        TypeError
            If the logfile is not of the expected type.

        Notes
        -----
        - The method modifies the log record in place.
        - It uses a dictionary `logmessage_to_formatter_map` to map specific log messages to formatter functions.
        - It suppresses log events listed in `suppress_events`.
        - It suppresses warning and info log messages based on the content of `suppress_warning_logmessages` and `supress_info_logmessages`.
        - If the log message contains "Complete log(s):", it overrides the log file.
        - It handles specific log events such as "run_info", "job_started", "job_info", and "job_error".
        """
        log_record = record.__dict__

        log_level = log_record.get("levelname", None)
        log_message = log_record.get("msg", None)
        log_event = log_record.get("event", None)
        execjob_logfile = log_record.get("log", None)

        if log_message is None and log_record.get("message") is not None:
            log_message = log_record.get("message")

        # log_record.get("log") will be either None or a list of the log files.
        # Since we always pass a single log file in the workflows we need to get the first element of the list.
        if execjob_logfile is not None:
            execjob_logfile = (
                execjob_logfile[0]
                if isinstance(execjob_logfile, list) and len(execjob_logfile) > 0
                else execjob_logfile
            )

        if log_level is None or log_message is None or log_message == "None":
            # If the record does not have a level or message, we cannot process it.
            return None

        # should lift this out to a better place to store this.
        # keeping it here for now to have an overview during development/migration.
        suppress_events = [
            "workflow_started",
            "resources_info",
            "shellcmd",
        ]
        suppress_warning_logmessages: list[str] = [
            "Your conda installation is not configured to use strict channel priorities.",
            "Failed to download resource needed for report:",
        ]
        supress_info_logmessages: list[str] = [
            "host:",
            "Assuming unrestricted shared filesystem usage.",
            "Using shell:",
            "Conda environments: ignored",
        ]
        logmessage_to_formatter_map: dict[str, Callable] = {
            "Activating conda environment": ColorizeLogMessagePath,
            "Activating singularity image": ColorizeLogMessagePath,
            "Creating conda environment": ColorizeLogMessagePath,
            "Removing incomplete Conda environment": ColorizeLogMessagePath,
            "Downloading and installing remote packages": CondaEnvInstallerPreamble,
            "Environment for": CondaEnvInstallSuccess,
            "Terminating process on user request": BaseSnakemakeAbortMessage,
            "Report created: results/snakemake_report.html.": ColorizeLogMessagePath,
            "Terminating processes on user request, this might take some time.": BaseSnakemakeAbortMessage,
            "will be created.": ColorizeLogMessagePath,
            # "Submitted DRMAA job": SubmitDRMAAMessage,
        }

        if log_event in suppress_events:
            return None

        if log_message is not None and any(
            x in log_message for x in suppress_warning_logmessages
        ):
            return None

        if log_message is not None and any(
            x in log_message for x in supress_info_logmessages
        ):
            return None

        if log_message is not None and any(
            x in log_message for x in list(logmessage_to_formatter_map.keys())
        ):
            for key in logmessage_to_formatter_map.keys():
                if key in log_message:
                    log_record["msg"] = logmessage_to_formatter_map[key](log_message)
                    return logging.makeLogRecord(log_record)

        if log_message is not None and "Complete log(s):" in log_message:
            log_record["msg"] = LogFileOverride(log_record, logfile)
            return logging.makeLogRecord(log_record)
        if log_event == "run_info":
            # print(log_message)
            log_record["msg"] = print_jobstatistics_logmessage(log_message)
            return logging.makeLogRecord(log_record)
        if log_event == "job_started":
            log_record = HandleJobStartedMessage(log_record)
            return logging.makeLogRecord(log_record)
        if log_event == "job_info":
            log_record = HandleJobInfoMessage(log_record)
            return logging.makeLogRecord(log_record)
        if log_event == "job_error":
            log_record = HandleJobErrorMessage(log_record)
            return logging.makeLogRecord(log_record)

        record = logging.makeLogRecord(log_record)
        return record


def CondaEnvInstallerPreamble(msg: str) -> str:
    logmessage = f"Creating conda environment: [yellow]{msg}[/yellow]"
    return logmessage


def CondaEnvInstallSuccess(_msg) -> str:
    return "[green]Conda environment created![/green]"


def BaseSnakemakeAbortMessage(msg: str) -> str:
    return f"[red]{msg}[/red]"


def ColorizeLogMessagePath(msg: str) -> str:
    """This function colorizes any file paths found in a log message with magenta color.

    Parameters
    ----------
    msg : dict[str, Any]
        The `msg` parameter is a dictionary that contains information about a log message.

    """
    msg = msg.strip(".").strip()
    messageparts = msg.split(" ")
    for i, part in enumerate(messageparts):
        part = part.replace("...", "") if "..." in part else part
        if pathlib.Path(part).exists():
            # If the part is a valid path, colorize it
            messageparts[i] = f"[magenta]{part}[/magenta]"
    # Join the parts back together into a single string
    formatted_message = " ".join(messageparts)
    return formatted_message


def LogFileOverride(msg, logfile) -> str:
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
    return f"Complete log: [magenta]{logfile}[/magenta]"


# TODO: update the DRMAA message handler to corrrectly return the properly formatted message.
# def SubmitDRMAAMessage(msg: dict[str, Any]) -> None:
#     """The function extracts job identifiers from a message dictionary and logs a message indicating that a
#     job has been submitted to a DRMAA cluster using an external job ID.

#     Parameters
#     ----------
#     msg : dict[str, Any]
#         The parameter `msg` is a dictionary that contains a log message.

#     """
#     if logmessage := msg.get("msg"):
#         logmessage_items = logmessage.rstrip(".").split(" ")
#         identifiers = []
#         for item in logmessage_items:
#             try:
#                 identifiers.append(int(item))
#             except Exception:
#                 continue
#         log.info(
#             f"Submitted JobID [cyan]{identifiers[0]}[/cyan] to the DRMAA cluster using external JobID [cyan]{identifiers[1]}[/cyan]."
#         )


def HandleJobStartedMessage(record: dict[str, Any]) -> dict[str, Any]:
    jobids: list = record.get("jobs", [])
    if len(jobids) == 1:
        record["msg"] = (
            f"Starting [cyan]{len(jobids)}[/cyan] task with jobid [cyan]{jobids[0]}[/cyan]"
        )
        return record
    # multiple jobs being started
    record["msg"] = (
        f"Starting [cyan]{len(jobids)}[/cyan] tasks with jobids: [cyan]{'[/cyan], [cyan]'.join(map(str, jobids))}[/cyan]"
    )
    return record


def HandleJobInfoMessage(record: dict[str, Any]) -> dict[str, Any]:
    """The function logs information about a job being executed, including the process name, sample, input
    target, reference ID, job ID, and log file.

    Parameters
    ----------
    msg : dict[str, Any]
        The `msg` parameter is a dictionary that contains information about a job to be executed.

    """
    process_name = record.get("rule_name", None)
    wildcards: dict[str, str] | None = record.get("wildcards", None)
    sample = wildcards.get("sample", None) if wildcards else None
    input_target = wildcards.get("Virus", None) if wildcards else None
    refid = wildcards.get("RefID", None) if wildcards else None
    list_logfile: list = record.get("log", [])
    job_id = record.get("jobid", None)
    logfile: str | None = (
        str(pathlib.Path(list_logfile[0]).absolute()) if list_logfile else None
    )
    if record.get("local", False):
        record["msg"] = (
            f"Executing localjob [green underline]{process_name}[/green underline] for sample [blue]{sample}[/blue] with target [blue]{input_target}[/blue] and reference-id [blue]{refid}[/blue]\nJob is using jobID [cyan]{job_id}[/cyan], logging output will be written to [magenta]{logfile}[/magenta]"
        )
        return record
    record["msg"] = (
        f"Executing job [green underline]{process_name}[/green underline] for sample [blue]{sample}[/blue] with target [blue]{input_target}[/blue] and reference-id [blue]{refid}[/blue]\nJob is using jobID [cyan]{job_id}[/cyan], logging output will be written to [magenta]{logfile}[/magenta]"
    )
    return record


def HandleJobErrorMessage(record: dict[str, Any]) -> dict[str, Any]:
    """This function handles error messages for a job and logs relevant information.

    Parameters
    ----------
    msg : dict[str, Any]
        The `msg` parameter is a dictionary containing information about a job error message.

    """
    process_name = record.get("rule_name", None)
    wildcards: dict[str, str] | None = record.get("wildcards", None)
    sample = wildcards.get("sample", None) if wildcards else None
    input_target = wildcards.get("Virus", None) if wildcards else None
    refid = wildcards.get("RefID", None) if wildcards else None
    list_logfile: list = record.get("log", [])
    job_id = record.get("jobid", None)
    logfile: str | None = (
        str(pathlib.Path(list_logfile[0]).absolute()) if list_logfile else None
    )
    shellcmd = record.get("shellcmd", "").strip().replace("\n", " ").replace("  ", " ")
    if outputfiles := record.get("output"):
        outputfiles_list = " ".join(list(outputfiles))
    else:
        outputfiles_list = "None"
    record["msg"] = (
        f"Job [red underline]{process_name}[/red underline] for sample [blue]{sample}[/blue] with target [blue]{input_target}[/blue] and reference-id [blue]{refid}[/blue] failed!\nJob is using jobID [cyan]{job_id}[/cyan], logging output will be written to [magenta]{logfile}[/magenta]\nThe expected output file(s) are as follows:\n[red]{outputfiles_list}[/red]\nThe following shell command was issued:\n[red]{shellcmd}[/red]"
    )
    return record


def print_jobstatistics_logmessage(msg: str) -> str:
    # if logmessage := msg.get("msg"):
    logmessage = msg.split("\n", 1)[1]
    return f"Workflow statistics:\n[yellow]{logmessage}[/yellow]"


log = logging.getLogger(__prog__)
log.propagate = False
log.setLevel(logging.INFO)
