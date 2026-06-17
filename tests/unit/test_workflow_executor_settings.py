from unittest import mock

import pytest

from ViroConstrictor.scheduler import Scheduler


class TestExecutorSettingsFor:

    @pytest.mark.parametrize(
        "scheduler",
        [Scheduler.LOCAL, Scheduler.DRYRUN, Scheduler.LSF],
    )
    def test_returns_none_for_in_process_or_unhandled(self, scheduler):
        from ViroConstrictor.workflow_executor import _executor_settings_for

        assert _executor_settings_for(scheduler) is None

    def test_returns_slurm_executor_settings(self):
        slurm_module = pytest.importorskip("snakemake_executor_plugin_slurm")
        from ViroConstrictor.workflow_executor import _executor_settings_for

        settings = _executor_settings_for(Scheduler.SLURM)
        assert isinstance(settings, slurm_module.ExecutorSettings)
        assert hasattr(settings, "logdir")
        assert hasattr(settings, "partition_config")


class TestWorkflowExecutorPassesSettings:
    """WorkflowExecutor must forward executor_settings into execute_workflow."""

    def _run_executor(self, scheduler):
        """Run WorkflowExecutor against a mocked SnakemakeApi, return the execute_workflow call args."""
        from ViroConstrictor.workflow_executor import WorkflowExecutor

        with mock.patch("ViroConstrictor.workflow_executor.SnakemakeApi") as api_cls:
            api_ctx = api_cls.return_value.__enter__.return_value
            dag_api = api_ctx.workflow.return_value.dag.return_value

            WorkflowExecutor(
                parsed_input=mock.MagicMock(),
                workflow_config=mock.MagicMock(),
                scheduler=scheduler,
            )

            return dag_api.execute_workflow.call_args

    def test_slurm_forwards_non_none_executor_settings(self):
        pytest.importorskip("snakemake_executor_plugin_slurm")
        call = self._run_executor(Scheduler.SLURM)
        assert call is not None, "execute_workflow was never called"
        assert call.kwargs.get("executor_settings") is not None

    def test_local_forwards_none_executor_settings(self):
        call = self._run_executor(Scheduler.LOCAL)
        assert call is not None
        # Absent or explicitly None are both fine here.
        assert call.kwargs.get("executor_settings") is None
