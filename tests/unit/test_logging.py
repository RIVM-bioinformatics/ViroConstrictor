"""Unit tests for ViroConstrictor logging formatters and handler dispatch."""

import contextlib
import datetime
import logging
from pathlib import Path
from typing import Any, Generator

import pytest

import ViroConstrictor.logging as vc_logging


@pytest.fixture
def restore_logger_handlers() -> Generator[None, None, None]:
    """Fixture to restore logger handlers and level after each test.

    Saves the current logger's handlers and level before the test, then
    restores them after the test completes to prevent state leakage between
    tests.

    Yields
    ------
    None
    """
    original_handlers = list(vc_logging.log.handlers)
    original_level = vc_logging.log.level
    try:
        yield
    finally:
        for handler in list(vc_logging.log.handlers):
            with contextlib.suppress(Exception):
                handler.close()
            vc_logging.log.removeHandler(handler)
        for handler in original_handlers:
            vc_logging.log.addHandler(handler)
        vc_logging.log.setLevel(original_level)


def test_setup_logger_creates_dir_and_handler(monkeypatch: pytest.MonkeyPatch, tmp_path: Path, restore_logger_handlers: None) -> None:
    """Test that setup_logger creates a directory and log handler.

    Verifies that setup_logger() creates the log directory if it doesn't
    exist and initializes a handler with a timestamped log file path.

    Parameters
    ----------
    monkeypatch : pytest.MonkeyPatch
        Pytest fixture for mocking attributes.
    tmp_path : Path
        Temporary directory provided by pytest.
    restore_logger_handlers : None
        Fixture to restore logger state after test.
    """
    captured: dict[str, str] = {}

    class _FakeDateTime:
        @staticmethod
        def now() -> datetime.datetime:
            return datetime.datetime(2024, 1, 2, 3, 4, 5)

    class _FakeDateTimeModule:
        datetime = _FakeDateTime

    def fake_handler_init(logfile_path: str, *args: Any, **kwargs: Any) -> None:
        captured["logfile"] = logfile_path

    monkeypatch.setattr(vc_logging, "datetime", _FakeDateTimeModule)
    monkeypatch.setattr(vc_logging, "ViroConstrictorBaseLogHandler", fake_handler_init)

    workdir = tmp_path / "new_logs"
    logfile = vc_logging.setup_logger(str(workdir))

    assert workdir.exists()
    assert logfile == captured["logfile"]
    assert "ViroConstrictor_2024-01-02T03:04:05.log" in logfile


@pytest.mark.parametrize(
    "level, expected_tabs",
    [
        ("DEBUG", "\n\t\t\t\t"),
        ("INFO", "\n\t\t\t\t\t\t\t"),
    ],
)
def test_strip_brackets_filter_removes_markup_and_formats_newlines(level: str, expected_tabs: str) -> None:
    """Test that StripBracketsFilter removes Rich markup and formats newlines.

    Verifies that StripBracketsFilter removes Rich markup tags like [red] and
    [/red], and replaces newlines with tabs appropriate to the log level.

    Parameters
    ----------
    level : str
        Log level name (DEBUG or INFO) - parametrized.
    expected_tabs : str
        Expected tab formatting for the log level - parametrized.
    """
    record = logging.makeLogRecord({"msg": "[red]hello[/red]\nworld", "levelname": level})

    filtered = vc_logging.StripBracketsFilter().filter(record)

    assert filtered is record
    assert "[red]" not in record.msg
    assert expected_tabs in record.msg


@pytest.mark.parametrize(
    "func,args,expected",
    [
        (vc_logging.CondaEnvInstallerPreamble, ("A",), "Creating conda environment: [yellow]A[/yellow]"),
        (vc_logging.CondaEnvInstallSuccess, ("ignored",), "[green]Conda environment created![/green]"),
        (vc_logging.BaseSnakemakeAbortMessage, ("Stop",), "[red]Stop[/red]"),
        (vc_logging.LogFileOverride, ({"msg": "Complete log(s):"}, "/tmp/x.log"), "Complete log: [magenta]/tmp/x.log[/magenta]"),
    ],
)
def test_simple_message_helpers(func: Any, args: tuple[Any, ...], expected: str) -> None:
    """Test message helper functions with various inputs.

    Verifies that simple message helper functions (CondaEnvInstallerPreamble,
    CondaEnvInstallSuccess, BaseSnakemakeAbortMessage, LogFileOverride)
    return correctly formatted Rich-markup strings.

    Parameters
    ----------
    func : Any
        Message helper function to test - parametrized.
    args : tuple[Any, ...]
        Arguments to pass to the function - parametrized.
    expected : str
        Expected string output - parametrized.
    """
    assert func(*args) == expected


def test_colorize_log_message_path_marks_existing_paths(tmp_path: Path) -> None:
    """Test that ColorizeLogMessagePath marks existing file paths.

    Verifies that ColorizeLogMessagePath() detects existing file paths in
    messages and wraps them with Rich magenta color markup.

    Parameters
    ----------
    tmp_path : Path
        Temporary directory provided by pytest.
    """
    existing = tmp_path / "exists.txt"
    existing.write_text("x", encoding="utf-8")

    message = f"Creating conda environment ... {existing}."
    result = vc_logging.ColorizeLogMessagePath(message)

    assert "[magenta]" in result
    assert str(existing) in result


def test_handle_job_started_message_single_and_multiple() -> None:
    """Test HandleJobStartedMessage formats singular and plural task counts.

    Verifies that HandleJobStartedMessage() correctly pluralizes "task" vs
    "tasks" depending on the number of started jobs.
    """
    single = vc_logging.HandleJobStartedMessage({"jobs": [7]})
    multiple = vc_logging.HandleJobStartedMessage({"jobs": [1, 2, 3]})

    assert "1[/cyan] task" in single["msg"]
    assert "3[/cyan] tasks" in multiple["msg"]


@pytest.mark.parametrize("local", [True, False])
def test_handle_job_info_message_variants(local: bool) -> None:
    """Test HandleJobInfoMessage formats job info with local/remote context.

    Verifies that HandleJobInfoMessage() correctly formats job information,
    including rule name, sample wildcards, job ID, and whether the job is
    executing locally or remotely.

    Parameters
    ----------
    local : bool
        Whether the job executes locally (True) or remotely (False) - parametrized.
    """
    record = {
        "rule_name": "align",
        "wildcards": {"sample": "S1", "Virus": "SARS", "RefID": "R1"},
        "log": ["/tmp/job.log"],
        "jobid": 11,
        "local": local,
    }

    out = vc_logging.HandleJobInfoMessage(record)

    assert "align" in out["msg"]
    assert "S1" in out["msg"]
    assert "jobID [cyan]11[/cyan]" in out["msg"]
    if local:
        assert "Executing localjob" in out["msg"]
    else:
        assert "Executing job" in out["msg"]


def test_handle_job_error_message_with_and_without_outputs() -> None:
    """Test HandleJobErrorMessage formats with and without job outputs.

    Verifies that HandleJobErrorMessage() correctly handles error records
    with and without output files, and normalizes shell command whitespace.
    """
    base_record = {
        "rule_name": "consensus",
        "wildcards": {"sample": "S1", "Virus": "SARS", "RefID": "R1"},
        "log": ["/tmp/err.log"],
        "jobid": 12,
        "shellcmd": " echo  hi\n there  ",
    }

    with_outputs = vc_logging.HandleJobErrorMessage({**base_record, "output": ["a.txt", "b.txt"]})
    without_outputs = vc_logging.HandleJobErrorMessage({**base_record, "output": []})

    assert "a.txt b.txt" in with_outputs["msg"]
    assert "None" in without_outputs["msg"]
    assert "echo hi there" in with_outputs["msg"]


def test_handle_job_debug_message_sets_debug_level() -> None:
    """Test that handle_job_debug_message sets log level to DEBUG.

    Verifies that handle_job_debug_message() converts a log record to DEBUG
    level and normalizes shell command whitespace.
    """
    record = {"msg": "orig", "levelname": "INFO", "levelno": 20}

    out = vc_logging.handle_job_debug_message(
        record=record,
        rule_name="r",
        wildcards={"sample": "S1"},
        rule_id="9",
        rule_input=["in1"],
        rule_output=["out1"],
        shellcmd="  cmd   --flag  ",
    )

    assert out["levelname"] == "DEBUG"
    assert out["levelno"] == 10
    assert "Snakemake workflow :: r" in out["msg"]
    assert "cmd --flag" in out["msg"]


def test_print_jobstatistics_logmessage() -> None:
    """Test that print_jobstatistics_logmessage formats workflow statistics.

    Verifies that print_jobstatistics_logmessage() wraps statistics content
    with a "Workflow statistics:" header and Rich yellow markup.
    """
    out = vc_logging.print_jobstatistics_logmessage("Header\nline1\nline2")
    assert out == "Workflow statistics:\n[yellow]line1\nline2[/yellow]"


def test_dispatch_returns_none_for_missing_required_fields(tmp_path: Path, restore_logger_handlers: None) -> None:
    """Test that dispatch returns None when required fields are missing.

    Verifies that _dispatch_log_record_formatter() returns None for records
    with missing required fields (e.g., None levelname).

    Parameters
    ----------
    tmp_path : Path
        Temporary directory provided by pytest.
    restore_logger_handlers : None
        Fixture to restore logger state after test.
    """
    handler = vc_logging.ViroConstrictorBaseLogHandler(str(tmp_path / "log.txt"))

    out = handler._dispatch_log_record_formatter(logging.makeLogRecord({"msg": "x", "levelname": None}))

    assert out is None


@pytest.mark.parametrize("event", ["workflow_started", "resources_info", "shellcmd"])
def test_dispatch_suppresses_events(event: str, tmp_path: Path, restore_logger_handlers: None) -> None:
    """Test that dispatch suppresses specific Snakemake events.

    Verifies that _dispatch_log_record_formatter() returns None for specific
    Snakemake events that should not be logged (workflow_started,
    resources_info, shellcmd).

    Parameters
    ----------
    event : str
        Event name to test - parametrized.
    tmp_path : Path
        Temporary directory provided by pytest.
    restore_logger_handlers : None
        Fixture to restore logger state after test.
    """
    handler = vc_logging.ViroConstrictorBaseLogHandler(str(tmp_path / "log.txt"))
    record = logging.makeLogRecord({"msg": "any", "levelname": "INFO", "event": event})

    assert handler._dispatch_log_record_formatter(record) is None


@pytest.mark.parametrize(
    "message",
    [
        "Your conda installation is not configured to use strict channel priorities.",
        "Failed to download resource needed for report:",
        "host:",
        "Assuming unrestricted shared filesystem usage.",
        "Using shell:",
        "Conda environments: ignored",
    ],
)
def test_dispatch_suppresses_specific_messages(message: str, tmp_path: Path, restore_logger_handlers: None) -> None:
    """Test that dispatch suppresses specific message patterns.

    Verifies that _dispatch_log_record_formatter() returns None for messages
    matching known suppression patterns (conda channel priorities, resource
    downloads, host info, etc.).

    Parameters
    ----------
    message : str
        Message text to test - parametrized.
    tmp_path : Path
        Temporary directory provided by pytest.
    restore_logger_handlers : None
        Fixture to restore logger state after test.
    """
    handler = vc_logging.ViroConstrictorBaseLogHandler(str(tmp_path / "log.txt"))
    record = logging.makeLogRecord({"msg": message, "levelname": "INFO"})

    assert handler._dispatch_log_record_formatter(record) is None


def test_dispatch_map_and_event_specific_branches(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, restore_logger_handlers: None) -> None:
    """Test dispatch mapping and event-specific formatter branches.

    Verifies that _dispatch_log_record_formatter() correctly handles:
    - Mapped messages (e.g., termination signals)
    - Special event types (run_info, job_started, job_info, job_error)
    - Complete log overrides

    Parameters
    ----------
    tmp_path : Path
        Temporary directory provided by pytest.
    monkeypatch : pytest.MonkeyPatch
        Pytest fixture for mocking attributes.
    restore_logger_handlers : None
        Fixture to restore logger state after test.
    """
    handler = vc_logging.ViroConstrictorBaseLogHandler(str(tmp_path / "log.txt"))
    monkeypatch.setattr(vc_logging, "logfile", "/tmp/override.log", raising=False)

    formatter_rec = logging.makeLogRecord(
        {
            "msg": "Terminating process on user request",
            "levelname": "INFO",
        }
    )
    out_formatter = handler._dispatch_log_record_formatter(formatter_rec)
    assert out_formatter is not None
    assert out_formatter[0]["msg"] == "[red]Terminating process on user request[/red]"

    complete_log_rec = logging.makeLogRecord(
        {
            "msg": "Complete log(s): location",
            "levelname": "INFO",
        }
    )
    out_complete_log = handler._dispatch_log_record_formatter(complete_log_rec)
    assert out_complete_log is not None
    assert "Complete log:" in out_complete_log[0]["msg"]

    run_info_rec = logging.makeLogRecord(
        {
            "msg": "Stats\nJobs: 1",
            "levelname": "INFO",
            "event": "run_info",
        }
    )
    out_run_info = handler._dispatch_log_record_formatter(run_info_rec)
    assert out_run_info is not None
    assert "Workflow statistics:" in out_run_info[0]["msg"]

    start_rec = logging.makeLogRecord(
        {
            "msg": "start",
            "levelname": "INFO",
            "event": "job_started",
            "jobs": [3],
        }
    )
    out_start = handler._dispatch_log_record_formatter(start_rec)
    assert out_start is not None
    assert "Starting" in out_start[0]["msg"]

    info_rec = logging.makeLogRecord(
        {
            "msg": "info",
            "levelname": "INFO",
            "event": "job_info",
            "rule_name": "r",
            "wildcards": {"sample": "s", "Virus": "v", "RefID": "id"},
            "log": ["/tmp/a.log"],
            "jobid": 1,
        }
    )
    out_info = handler._dispatch_log_record_formatter(info_rec)
    assert out_info is not None
    assert "Executing job" in out_info[0]["msg"]

    err_rec = logging.makeLogRecord(
        {
            "msg": "err",
            "levelname": "INFO",
            "event": "job_error",
            "rule_name": "r",
            "wildcards": {"sample": "s", "Virus": "v", "RefID": "id"},
            "log": ["/tmp/a.log"],
            "jobid": 1,
            "shellcmd": "echo x",
            "output": ["x.txt"],
        }
    )
    out_err = handler._dispatch_log_record_formatter(err_rec)
    assert out_err is not None
    assert "failed!" in out_err[0]["msg"]


def test_dispatch_adds_debug_record_when_job_fields_are_present(tmp_path: Path, restore_logger_handlers: None) -> None:
    """Test that dispatch adds a DEBUG record when job-related fields exist.

    Verifies that _dispatch_log_record_formatter() adds a DEBUG-level record
    with job details (rule, sample, inputs, outputs, command) when these
    fields are present in the input record.

    Parameters
    ----------
    tmp_path : Path
        Temporary directory provided by pytest.
    restore_logger_handlers : None
        Fixture to restore logger state after test.
    """
    handler = vc_logging.ViroConstrictorBaseLogHandler(str(tmp_path / "log.txt"))

    record = logging.makeLogRecord(
        {
            "msg": "plain",
            "levelname": "INFO",
            "rule_name": "r",
            "wildcards": {"sample": "S1"},
            "jobid": 2,
            "input": ["i"],
            "output": ["o"],
            "shellcmd": "cmd",
        }
    )
    out = handler._dispatch_log_record_formatter(record)

    assert out is not None
    assert any(item["levelname"] == "DEBUG" for item in out)
    assert any(item["msg"] == "plain" for item in out)


def test_emit_filters_by_level_and_handles_errors(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, restore_logger_handlers: None) -> None:
    """Test that emit filters records by level and handles handler errors.

    Verifies that the handler's emit() method:
    - Filters records below the configured log level
    - Dispatches records to appropriate handlers
    - Catches and handles emitter exceptions without propagating

    Parameters
    ----------
    tmp_path : Path
        Temporary directory provided by pytest.
    monkeypatch : pytest.MonkeyPatch
        Pytest fixture for mocking attributes.
    restore_logger_handlers : None
        Fixture to restore logger state after test.
    """
    handler = vc_logging.ViroConstrictorBaseLogHandler(str(tmp_path / "log.txt"))

    emitted_console: list[logging.LogRecord] = []
    emitted_file: list[logging.LogRecord] = []
    handled_errors: list[logging.LogRecord] = []

    def fake_console_emit(record: logging.LogRecord) -> None:
        if "raise-console" in str(record.msg):
            raise RuntimeError("console")
        emitted_console.append(record)

    def fake_file_emit(record: logging.LogRecord) -> None:
        if "raise-file" in str(record.msg):
            raise RuntimeError("file")
        emitted_file.append(record)

    monkeypatch.setattr(handler.console_handler, "emit", fake_console_emit)
    monkeypatch.setattr(handler.instance_file_handler, "emit", fake_file_emit)
    monkeypatch.setattr(handler, "handleError", lambda rec: handled_errors.append(rec))

    monkeypatch.setattr(
        handler,
        "_dispatch_log_record_formatter",
        lambda _rec: [
            {"msg": "below-level", "levelname": "DEBUG", "levelno": 10},
            {"msg": "ok", "levelname": "INFO", "levelno": 20},
            {"msg": "raise-console", "levelname": "INFO", "levelno": 20},
            {"msg": "raise-file", "levelname": "INFO", "levelno": 20},
        ],
    )

    handler.emit(logging.makeLogRecord({"msg": "base", "levelname": "INFO", "levelno": 20}))

    assert any(str(r.msg) == "ok" for r in emitted_console)
    assert all(str(r.msg) != "below-level" for r in emitted_console)
    assert any(str(r.msg) == "ok" for r in emitted_file)
    assert len(handled_errors) == 2


def test_emit_returns_when_dispatch_is_none(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, restore_logger_handlers: None) -> None:
    """Test that emit returns early when dispatch returns None.

    Verifies that when _dispatch_log_record_formatter() returns None, the
    emit() method does not dispatch to any handlers.

    Parameters
    ----------
    tmp_path : Path
        Temporary directory provided by pytest.
    monkeypatch : pytest.MonkeyPatch
        Pytest fixture for mocking attributes.
    restore_logger_handlers : None
        Fixture to restore logger state after test.
    """
    handler = vc_logging.ViroConstrictorBaseLogHandler(str(tmp_path / "log.txt"))
    called: dict[str, bool] = {"console": False, "file": False}

    monkeypatch.setattr(handler, "_dispatch_log_record_formatter", lambda _rec: None)
    monkeypatch.setattr(handler.console_handler, "emit", lambda _rec: called.__setitem__("console", True))
    monkeypatch.setattr(handler.instance_file_handler, "emit", lambda _rec: called.__setitem__("file", True))

    handler.emit(logging.makeLogRecord({"msg": "ignored", "levelname": "INFO", "levelno": 20}))

    assert called == {"console": False, "file": False}


def test_close_closes_handlers(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, restore_logger_handlers: None) -> None:
    """Test that close() method closes both console and file handlers.

    Verifies that the handler's close() method properly delegates to and
    closes both the console handler and file handler.

    Parameters
    ----------
    tmp_path : Path
        Temporary directory provided by pytest.
    monkeypatch : pytest.MonkeyPatch
        Pytest fixture for mocking attributes.
    restore_logger_handlers : None
        Fixture to restore logger state after test.
    """
    handler = vc_logging.ViroConstrictorBaseLogHandler(str(tmp_path / "log.txt"))
    closed = {"console": False, "file": False}

    monkeypatch.setattr(handler.console_handler, "close", lambda: closed.__setitem__("console", True))
    monkeypatch.setattr(handler.instance_file_handler, "close", lambda: closed.__setitem__("file", True))

    handler.close()

    assert closed == {"console": True, "file": True}
