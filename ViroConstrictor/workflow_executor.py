from pathlib import Path

from snakemake.api import SnakemakeApi
from snakemake.logging import logger, logger_manager
from snakemake_interface_common.exceptions import WorkflowError

from ViroConstrictor.logging import log
from ViroConstrictor.parser import CLIparser
from ViroConstrictor.scheduler import Scheduler
from ViroConstrictor.workflow_config import WorkflowConfig


# TODO: Modify this method to re-add remote execution settings if they are set in the CLIparser or config file.
def run_snakemake_workflow(inputs_obj: CLIparser, stage: str, scheduler: Scheduler) -> tuple[bool, WorkflowConfig]:
    """
    Run the snakemake workflow for the specified stage.

    Parameters
    ----------
    inputs_obj : CLIparser
        The parsed command line arguments.
    stage : str
        The stage of the workflow to run.

    Returns
    -------
    bool
        True if the workflow ran successfully, False otherwise.
    """
    workflow_configuration = WorkflowConfig(parsed_inputs=inputs_obj, outdir_override="", vc_stage=stage)
    try:
        WorkflowExecutor(
            parsed_input=inputs_obj,
            workflow_config=workflow_configuration,
            scheduler=scheduler,
        )
    except WorkflowError as e:
        log.error(f"Workflow execution failed with error: {e}\nPlease check the logs and your settings for more information and try again later.")
        return False, workflow_configuration
    return True, workflow_configuration


class WorkflowExecutor:
    """
    Executes a Snakemake workflow based on the provided configuration.
    Attributes
    ----------
    parsed_input : CLIparser
        The parsed command-line input arguments.
    workflow_config : WorkflowConfig
        The configuration settings for the workflow.
    workflow_api : SnakemakeWorkflowApi
        The Snakemake workflow API instance.
    dag_api : SnakemakeDagApi
        The Snakemake DAG API instance.
    Methods
    -------
    __init__(parsed_input: CLIparser, workflow_config: WorkflowConfig)
        Initializes the WorkflowExecutor, configures the Snakemake workflow and DAG,
        and executes the workflow.
    """

    def __init__(
        self,
        parsed_input: CLIparser,
        workflow_config: WorkflowConfig,
        scheduler: Scheduler = Scheduler.LOCAL,
    ) -> None:
        self.parsed_input = parsed_input
        self.workflow_config = workflow_config

        with SnakemakeApi(output_settings=self.workflow_config.output_settings) as snakemake_api:
            self.workflow_api = snakemake_api.workflow(
                resource_settings=self.workflow_config.resource_settings,
                config_settings=self.workflow_config.workflow_configsettings,
                storage_settings=self.workflow_config.storage_settings,
                workflow_settings=self.workflow_config.workflow_settings,
                deployment_settings=self.workflow_config.deployment_settings,
                snakefile=Path(workflow_config.workflow_file),
                workdir=Path(parsed_input.workdir),
            )

            self.dag_api = self.workflow_api.dag(self.workflow_config.dag_settings)

            self.dag_api.execute_workflow(
                executor=scheduler.value[0],
                execution_settings=self.workflow_config.execution_settings,
                remote_execution_settings=self.workflow_config.remote_execution_settings,
                scheduling_settings=self.workflow_config.scheduling_settings,
            )

        # 'Forcefully' stop snakemake's logger manager
        # This is necessary as otherwise the logger_manager will stay initialized
        # which will block it from being initialized again in a subsequent workflow call (i.e. first match_ref, then main)
        # and because the logger_manager will not be re-initialized, our custom logging plugin will not be used.
        if logger_manager.queue_listener is not None:
            logger_manager.stop()

        # Reset the initialized flag to allow the setup() method to run again
        logger_manager.initialized = False

        # Remove the QueueHandler that was added during setup() to avoid
        # duplicate handlers on the next initialization.
        for handler in list(logger.handlers):
            logger.removeHandler(handler)
