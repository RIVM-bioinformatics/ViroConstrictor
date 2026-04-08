"""Test workflow execution setup and Snakemake API integration.

These tests cover debugger patching, workflow construction, execution error
handling, and logger cleanup around the executor lifecycle.
"""

from pathlib import Path
from types import SimpleNamespace
from typing import Any, cast

import pytest
from snakemake_interface_common.exceptions import WorkflowError

import ViroConstrictor.workflow_executor as workflow_executor
from ViroConstrictor.parser import CLIparser
from ViroConstrictor.scheduler import Scheduler
from ViroConstrictor.workflow_config import WorkflowConfig


class _FakeSnakemakeApiContext:
    """Wrap a fake Snakemake API in a context manager."""

    def __init__(self, snakemake_api: Any):
        self._snakemake_api = snakemake_api

    def __enter__(self) -> Any:
        return self._snakemake_api

    def __exit__(self, exc_type: Any, exc_value: Any, traceback: Any) -> bool:
        return False


class _FakeLogger:
    """Provide a lightweight logger stand-in for handler cleanup tests."""

    def __init__(self, handlers: list[Any] | None = None) -> None:
        self.handlers = handlers if handlers is not None else []
        self.removed_handlers: list[Any] = []

    def removeHandler(self, handler: Any) -> None:  # noqa: N802 - external API name
        self.removed_handlers.append(handler)
        self.handlers.remove(handler)


class _FakeLoggerManager:
    """Track logger-manager stop calls during executor teardown."""

    def __init__(self, queue_listener: Any) -> None:
        self.queue_listener = queue_listener
        self.initialized = True
        self.stop_call_count = 0

    def stop(self) -> None:
        self.stop_call_count += 1


def _make_parsed_input(tmp_path: Path) -> SimpleNamespace:
    """Build parsed input state with a temporary workdir.

    Parameters
    ----------
    tmp_path : Path
        Temporary directory used to derive the workflow working directory.

    Returns
    -------
    SimpleNamespace
        Parsed input object exposing the workdir attribute used by the executor.
    """
    return SimpleNamespace(workdir=str(tmp_path / "wd"))


def _make_workflow_config() -> SimpleNamespace:
    """Build a workflow config stand-in with executor-facing attributes."""
    return SimpleNamespace(
        output_settings=object(),
        resource_settings=object(),
        workflow_configsettings=object(),
        storage_settings=object(),
        workflow_settings=object(),
        deployment_settings=object(),
        workflow_file="main/workflow.smk",
        dag_settings=object(),
        execution_settings=object(),
        remote_execution_settings=object(),
        scheduling_settings=object(),
    )


def test_patch_debugger_returns_early_without_active_debugger(monkeypatch: pytest.MonkeyPatch) -> None:
    """Verify _patch_debugger_for_snakemake returns early when no debugger is active."""
    monkeypatch.setattr(workflow_executor.sys, "gettrace", lambda: None)

    def original_abs_and_canonical(filename: Any, norm_paths_container: Any) -> tuple[Any, Any]:
        return filename, norm_paths_container

    fake_utils = SimpleNamespace(
        _abs_and_canonical_path=original_abs_and_canonical,
        NORM_PATHS_CONTAINER={"k": "v"},
    )
    monkeypatch.setitem(workflow_executor.sys.modules, "pydevd_file_utils", fake_utils)

    workflow_executor._patch_debugger_for_snakemake()

    assert fake_utils._abs_and_canonical_path is original_abs_and_canonical


def test_patch_debugger_returns_when_pydevd_module_missing(monkeypatch: pytest.MonkeyPatch) -> None:
    """Verify _patch_debugger_for_snakemake handles missing pydevd_file_utils module gracefully."""
    monkeypatch.setattr(workflow_executor.sys, "gettrace", lambda: object())
    monkeypatch.delitem(workflow_executor.sys.modules, "pydevd_file_utils", raising=False)

    workflow_executor._patch_debugger_for_snakemake()


def test_patch_debugger_wraps_annotated_string_filename(monkeypatch: pytest.MonkeyPatch) -> None:
    """Verify _patch_debugger_for_snakemake correctly wraps annotated string filename handling."""
    monkeypatch.setattr(workflow_executor.sys, "gettrace", lambda: object())

    captured: dict[str, Any] = {}

    def original_abs_and_canonical(filename: Any, norm_paths_container: Any) -> tuple[Any, Any]:
        captured["filename"] = filename
        captured["norm"] = norm_paths_container
        return filename, norm_paths_container

    fake_utils = SimpleNamespace(
        _abs_and_canonical_path=original_abs_and_canonical,
        NORM_PATHS_CONTAINER={"k": "v"},
    )
    monkeypatch.setitem(workflow_executor.sys.modules, "pydevd_file_utils", fake_utils)

    workflow_executor._patch_debugger_for_snakemake()

    class AnnotatedString:
        def __str__(self) -> str:
            return "/tmp/workflow.smk"

    result = fake_utils._abs_and_canonical_path(AnnotatedString())

    assert captured["filename"] == "/tmp/workflow.smk"
    assert captured["norm"] == {"k": "v"}
    assert result == ("/tmp/workflow.smk", {"k": "v"})


def test_run_snakemake_workflow_happy_path(monkeypatch: pytest.MonkeyPatch) -> None:
    """Verify run_snakemake_workflow initializes config and executor with correct arguments."""
    inputs = SimpleNamespace()
    scheduler = Scheduler.SLURM

    expected_config = SimpleNamespace(name="cfg")
    calls: dict[str, Any] = {}

    def fake_workflow_config(parsed_inputs: Any, outdir_override: str, vc_stage: str) -> Any:
        calls["config_args"] = (parsed_inputs, outdir_override, vc_stage)
        return expected_config

    def fake_executor(parsed_input: Any, workflow_config: Any, scheduler: Any) -> None:
        calls["executor_args"] = (parsed_input, workflow_config, scheduler)

    monkeypatch.setattr(workflow_executor, "WorkflowConfig", fake_workflow_config)
    monkeypatch.setattr(workflow_executor, "WorkflowExecutor", fake_executor)

    ok, cfg = workflow_executor.run_snakemake_workflow(inputs_obj=cast(CLIparser, inputs), stage="MAIN", scheduler=scheduler)

    assert ok is True
    assert cfg is expected_config
    assert calls["config_args"] == (inputs, "", "MAIN")
    assert calls["executor_args"] == (inputs, expected_config, scheduler)


def test_run_snakemake_workflow_handles_workflow_error(monkeypatch: pytest.MonkeyPatch) -> None:
    """Verify run_snakemake_workflow catches and reports WorkflowError from executor."""
    inputs = SimpleNamespace()
    expected_config = SimpleNamespace(name="cfg")

    def fake_workflow_config(parsed_inputs: Any, outdir_override: str, vc_stage: str) -> Any:
        return expected_config

    def fake_executor(parsed_input: Any, workflow_config: Any, scheduler: Any) -> None:
        raise WorkflowError("boom")

    captured_errors: list[str] = []
    monkeypatch.setattr(workflow_executor, "WorkflowConfig", fake_workflow_config)
    monkeypatch.setattr(workflow_executor, "WorkflowExecutor", fake_executor)
    monkeypatch.setattr(workflow_executor.log, "error", lambda msg: captured_errors.append(msg))

    ok, cfg = workflow_executor.run_snakemake_workflow(inputs_obj=cast(CLIparser, inputs), stage="MR", scheduler=Scheduler.LOCAL)

    assert ok is False
    assert cfg is expected_config
    assert len(captured_errors) == 1
    assert "Workflow execution failed with error: boom" in captured_errors[0]
    assert "Please check the logs and your settings" in captured_errors[0]


def test_workflow_executor_executes_workflow_and_cleans_logger_state(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """Verify WorkflowExecutor executes workflow and cleans up logger handlers on completion."""
    parsed_input = _make_parsed_input(tmp_path)
    workflow_config = _make_workflow_config()

    dag_calls: dict[str, Any] = {}
    workflow_calls: dict[str, Any] = {}

    class FakeDagApi:
        def execute_workflow(self, **kwargs: Any) -> None:
            dag_calls["execute_kwargs"] = kwargs

    class FakeWorkflowApi:
        def dag(self, dag_settings: Any) -> FakeDagApi:
            workflow_calls["dag_settings"] = dag_settings
            return FakeDagApi()

    class FakeSnakemakeApi:
        def workflow(self, **kwargs: Any) -> FakeWorkflowApi:
            workflow_calls["workflow_kwargs"] = kwargs
            return FakeWorkflowApi()

    snakemake_constructor_calls: dict[str, Any] = {}

    def fake_snakemake_api_constructor(output_settings: Any) -> _FakeSnakemakeApiContext:
        snakemake_constructor_calls["output_settings"] = output_settings
        return _FakeSnakemakeApiContext(FakeSnakemakeApi())

    fake_logger = _FakeLogger(handlers=[object(), object()])
    fake_logger_manager = _FakeLoggerManager(queue_listener=object())

    monkeypatch.setattr(workflow_executor, "SnakemakeApi", fake_snakemake_api_constructor)
    monkeypatch.setattr(workflow_executor, "logger", fake_logger)
    monkeypatch.setattr(workflow_executor, "logger_manager", fake_logger_manager)

    executor = workflow_executor.WorkflowExecutor(
        parsed_input=cast(CLIparser, parsed_input),
        workflow_config=cast(WorkflowConfig, workflow_config),
        scheduler=Scheduler.SLURM,
    )

    assert executor.parsed_input is parsed_input
    assert executor.workflow_config is workflow_config

    assert snakemake_constructor_calls["output_settings"] is workflow_config.output_settings

    assert workflow_calls["workflow_kwargs"]["resource_settings"] is workflow_config.resource_settings
    assert workflow_calls["workflow_kwargs"]["config_settings"] is workflow_config.workflow_configsettings
    assert workflow_calls["workflow_kwargs"]["storage_settings"] is workflow_config.storage_settings
    assert workflow_calls["workflow_kwargs"]["workflow_settings"] is workflow_config.workflow_settings
    assert workflow_calls["workflow_kwargs"]["deployment_settings"] is workflow_config.deployment_settings
    assert workflow_calls["workflow_kwargs"]["snakefile"] == Path("main/workflow.smk")
    assert workflow_calls["workflow_kwargs"]["workdir"] == Path(parsed_input.workdir)

    assert workflow_calls["dag_settings"] is workflow_config.dag_settings

    assert dag_calls["execute_kwargs"]["executor"] == "slurm"
    assert dag_calls["execute_kwargs"]["execution_settings"] is workflow_config.execution_settings
    assert dag_calls["execute_kwargs"]["remote_execution_settings"] is workflow_config.remote_execution_settings
    assert dag_calls["execute_kwargs"]["scheduling_settings"] is workflow_config.scheduling_settings

    assert fake_logger_manager.stop_call_count == 1
    assert fake_logger_manager.initialized is False
    assert fake_logger.handlers == []
    assert len(fake_logger.removed_handlers) == 2


def test_workflow_executor_does_not_stop_logger_when_queue_listener_absent(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """Verify WorkflowExecutor skips logger cleanup when queue_listener is not initialized."""
    parsed_input = _make_parsed_input(tmp_path)
    workflow_config = _make_workflow_config()

    class FakeDagApi:
        def execute_workflow(self, **kwargs: Any) -> None:
            return None

    class FakeWorkflowApi:
        def dag(self, dag_settings: Any) -> FakeDagApi:
            return FakeDagApi()

    class FakeSnakemakeApi:
        def workflow(self, **kwargs: Any) -> FakeWorkflowApi:
            return FakeWorkflowApi()

    fake_logger = _FakeLogger(handlers=[])
    fake_logger_manager = _FakeLoggerManager(queue_listener=None)

    monkeypatch.setattr(workflow_executor, "SnakemakeApi", lambda output_settings: _FakeSnakemakeApiContext(FakeSnakemakeApi()))
    monkeypatch.setattr(workflow_executor, "logger", fake_logger)
    monkeypatch.setattr(workflow_executor, "logger_manager", fake_logger_manager)

    workflow_executor.WorkflowExecutor(
        parsed_input=cast(CLIparser, parsed_input),
        workflow_config=cast(WorkflowConfig, workflow_config),
        scheduler=Scheduler.LOCAL,
    )

    assert fake_logger_manager.stop_call_count == 0
    assert fake_logger_manager.initialized is False


@pytest.mark.xfail(reason="Logger cleanup should still run when execute_workflow raises, but current implementation skips it")
def test_workflow_executor_cleans_logger_state_even_on_execute_failure(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """Verify WorkflowExecutor cleans up logger state even when DAG execution raises an exception."""
    parsed_input = _make_parsed_input(tmp_path)
    workflow_config = _make_workflow_config()

    class FakeDagApi:
        def execute_workflow(self, **kwargs: Any) -> None:
            raise RuntimeError("dag failed")

    class FakeWorkflowApi:
        def dag(self, dag_settings: Any) -> FakeDagApi:
            return FakeDagApi()

    class FakeSnakemakeApi:
        def workflow(self, **kwargs: Any) -> FakeWorkflowApi:
            return FakeWorkflowApi()

    fake_logger = _FakeLogger(handlers=[object()])
    fake_logger_manager = _FakeLoggerManager(queue_listener=object())

    monkeypatch.setattr(workflow_executor, "SnakemakeApi", lambda output_settings: _FakeSnakemakeApiContext(FakeSnakemakeApi()))
    monkeypatch.setattr(workflow_executor, "logger", fake_logger)
    monkeypatch.setattr(workflow_executor, "logger_manager", fake_logger_manager)

    with pytest.raises(RuntimeError, match="dag failed"):
        workflow_executor.WorkflowExecutor(
            parsed_input=cast(CLIparser, parsed_input),
            workflow_config=cast(WorkflowConfig, workflow_config),
            scheduler=Scheduler.LOCAL,
        )

    assert fake_logger_manager.stop_call_count == 1
    assert fake_logger_manager.initialized is False
    assert fake_logger.handlers == []
