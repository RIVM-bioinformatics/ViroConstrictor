import multiprocessing
import os
import sys
from argparse import Namespace
from configparser import ConfigParser
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Hashable

from snakemake.api import (
    ConfigSettings,
    DAGSettings,
    DeploymentMethod,
    DeploymentSettings,
    ExecutionSettings,
    OutputSettings,
    RemoteExecutionSettings,
    ResourceSettings,
    SchedulingSettings,
    StorageSettings,
    WorkflowSettings,
)
from snakemake.resources import DefaultResources
from snakemake.settings.enums import Quietness
from snakemake_interface_executor_plugins.settings import ExecMode
from snakemake_interface_logger_plugins.settings import (
    LogHandlerSettingsBase,
)

from ViroConstrictor.logging import log
from ViroConstrictor.parser import CLIparser
from ViroConstrictor.runconfigs import WriteYaml
from ViroConstrictor.workflow.containers import (
    construct_container_bind_args,
    download_containers,
)


def correct_unidirectional_flag(
    samples_dict: dict[Hashable, Any], flags: Namespace
) -> bool:
    """Corrects the unidirectional flag based on the platform and samples dictionary."""
    if flags.platform == "illumina" and flags.unidirectional is True:
        # check if the samples_dict has the INPUTFILE key, if it does, set the unidirectional flag to True
        if any("INPUTFILE" in sample for sample in samples_dict.values()):
            return True
        # check if the samples_dict has both R1 and R2 keys, if it does, set the unidirectional flag to False
        elif any("R1" in sample and "R2" in sample for sample in samples_dict.values()):
            log.warning(
                "Unidirectional flag is set to True, but both illumina R1 and R2 inputfiles were found. "
                "This may lead to unexpected behavior. Setting unidirectional flag to False."
            )
            return False
    return True  # Default to True if no specific conditions are met


class WorkflowConfig:
    """
    Configuration class for the ViroConstrictor workflow.

    Attributes
    ----------
    inputs : CLIparser
        Parsed command-line inputs.
    outdir_override : str
        Override for the output directory.
    configuration : dict
        User-defined configuration settings.
    vc_stage : str
        ViroConstrictor stage ('MR' or 'MAIN').
    samplesheetfilename : str
        Name of the samplesheet file.
    workflow_file : str
        Path to the Snakemake workflow file.
    output_settings : OutputSettings
        Settings related to output behavior.
    resource_settings : ResourceSettings
        Settings related to resource allocation.
    storage_settings : StorageSettings
        Settings related to storage.
    deployment_settings : DeploymentSettings
        Settings related to deployment (e.g., Conda, Apptainer).
    execution_settings : ExecutionSettings
        Settings related to Snakemake execution.
    scheduling_settings : SchedulingSettings
        Settings related to job scheduling.
    workflow_settings : WorkflowSettings
        Settings related to the overall workflow.
    snakemake_base_params : dict
        Base parameters passed to Snakemake.
    workflow_configsettings : ConfigSettings
        Settings related to workflow configuration.
    remote_execution_settings : RemoteExecutionSettings
        Settings for remote job execution.
    dag_settings : DAGSettings
        Settings related to DAG generation.

    Methods
    ----------
    _set_cores(cores: int) -> int
        Sets the number of cores to use, ensuring it doesn't exceed available cores.
    _get_max_local_mem() -> int
        Gets the maximum local memory available in MB.
    """

    def __init__(
        self,
        parsed_inputs: CLIparser,
        outdir_override: str = "",
        vc_stage: str | None = None,
    ) -> None:
        self.inputs = parsed_inputs
        self.outdir_override = outdir_override
        self.configuration = self.inputs.user_config
        self.vc_stage = vc_stage
        self.dryrun = self.inputs.flags.dryrun

        # check if VC_stage is set to either "MR" or "MAIN", these are the only two valid stages. Exit if None
        if self.vc_stage not in ["MR", "MAIN"]:
            raise ValueError("VC_stage must be set to either 'MR' or 'MAIN'.")
        if self.vc_stage == "MR":
            self.samplesheetfilename = "samples_mr"
            self.workflow_file = parsed_inputs.match_ref_snakefile
        if self.vc_stage == "MAIN":
            self.samplesheetfilename = "samples_main"
            self.workflow_file = parsed_inputs.snakefile

        assign_threads = MaxThreadsPerType(self.inputs, self.configuration)

        # TODO: dryrun only seems to work if outputsettings.dryrun is set to True, the executor is set to "dryrun" and the quietness is set to "SUBPROCESS"
        # Quite convoluted, but this is the only combination of settings that seems to actually run the workflow in true dryrun mode.
        self.output_settings = OutputSettings(
            dryrun=self.dryrun,  # this doesn't seem to actually work all that much? if dryrun is set to True, it will still run the workflow and actually execute the tasks.
            printshellcmds=False,
            nocolor=False,
            debug_dag=False,
            verbose=False,
            show_failed_logs=False,
            log_handler_settings={"viroconstrictor": LogHandlerSettingsBase()},
            keep_logger=False,
            stdout=False,
            benchmark_extended=False,
            quiet={Quietness.ALL},
        )

        # TODO: set resources dynamically based on the user configuration and the possible grid flags (remote executor plugins, like LSF, SLURM or other)
        # NOTE: current resource settings are just an example for LSF but this should be lifted to a separate function
        default_resource_setting = DefaultResources()
        default_resource_setting.set_resource(
            "lsf_queue", "bio"
        )  # this should be set dynamically based on the user config
        # default_resource_setting.set_resource("lsf_project", "PROJECT HERE") # this should be set dynamically as well, or be left empty if no default project can be found.

        self.resource_settings = ResourceSettings(
            cores=(
                300
                if self.configuration["COMPUTING"]["compmode"] == "grid"
                else self._set_cores(self.inputs.flags.threads)
            ),
            resources={"max_local_mem": self._get_max_local_mem()},
            nodes=200 if self.configuration["COMPUTING"]["compmode"] == "grid" else 1,
            default_resources=default_resource_setting,
        )

        self.storage_settings = StorageSettings()

        self.deployment_settings = DeploymentSettings(
            deployment_method=(
                {DeploymentMethod.APPTAINER}
                if self.configuration["REPRODUCTION"]["repro_method"] == "containers"
                else {DeploymentMethod.CONDA}
            ),
            conda_prefix=None,
            conda_cleanup_pkgs=None,
            conda_base_path=None,
            apptainer_prefix=(
                Path(self.configuration["REPRODUCTION"]["container_cache_path"])
                if self.configuration["REPRODUCTION"]["repro_method"] == "containers"
                else None
            ),
            apptainer_args=construct_container_bind_args(self.inputs.samples_dict),
        )

        self.execution_settings = ExecutionSettings(
            latency_wait=60,
            keep_going=True,
            debug=False,
            standalone=False,
            ignore_ambiguity=False,
            lock=True,
            ignore_incomplete=False,
            attempt=1,
            retries=3,
            use_threads=False,
        )

        self.scheduling_settings = SchedulingSettings(
            scheduler="greedy",
        )

        self.workflow_settings = WorkflowSettings(exec_mode=ExecMode.SUBPROCESS if self.dryrun else ExecMode.DEFAULT)

        unidirectional = correct_unidirectional_flag(
            self.inputs.samples_dict, self.inputs.flags
        )

        self.snakemake_base_params = {
            "sample_sheet": WriteYaml(
                self.inputs.samples_dict,
                f"{self.inputs.workdir}/{self.samplesheetfilename}.yaml",
            ),
            "logfile": self.inputs.logfile,
            "computing_execution": self.configuration["COMPUTING"]["compmode"],
            "max_local_mem": self._get_max_local_mem(),
            "platform": self.inputs.flags.platform,
            "unidirectional": unidirectional,
            "amplicon_type": self.inputs.flags.amplicon_type,
            "outdirOverride": self.outdir_override,
            "threads": {
                "Alignments": assign_threads.highcpu,
                "QC": assign_threads.midcpu,
                "AdapterRemoval": assign_threads.lowcpu,
                "PrimerRemoval": assign_threads.highcpu,
                "Consensus": assign_threads.midcpu,
                "Index": assign_threads.lowcpu,
                "Typing": assign_threads.lowcpu,
            },
        }

        self.workflow_configsettings = ConfigSettings(
            configfiles=[
                Path(
                    WriteYaml(
                        self.snakemake_base_params, f"{self.inputs.workdir}/config.yaml"
                    )
                )
            ],
        )

        self.remote_execution_settings = RemoteExecutionSettings(
            jobname="ViroConstrictor_{name}.jobid{jobid}",
            immediate_submit=False,
            envvars=[],
            max_status_checks_per_second=1.0,
        )

        self.dag_settings = DAGSettings(
            force_incomplete=True,
        )

        if self.deployment_settings.apptainer_prefix is not None:
            if (
                download_containers(
                    self.deployment_settings.apptainer_prefix,
                    self.output_settings.dryrun,
                )
                != 0
            ):
                log.error(
                    "Failed to download containers required for workflow.\nPlease check the logs and your settings for more information and try again later."
                )
                sys.exit(1)

    def _set_cores(self, cores: int) -> int:
        available: int = multiprocessing.cpu_count()
        if cores == available:
            return cores - 2
        return available - 2 if cores > available else cores

    def _get_max_local_mem(self) -> int:
        """Get the maximum local memory available in MB, minus a buffer of 2000 MB.

        Returns
        -------
        int
            The amount of memory available in MB.

        """
        avl_mem_bytes = os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES")
        return int(round(avl_mem_bytes / (1024.0**2) - 2000, -3))


@dataclass
class MaxThreadsPerType:
    """
    Represents the maximum number of threads to use for each type of process.

    Attributes
    ----------
    highcpu : int
        Maximum number of threads for high CPU processes.
    midcpu : int
        Maximum number of threads for medium CPU processes.
    lowcpu : int
        Maximum number of threads for low CPU processes.
    assignment : bool
        A flag to indicate if the computing mode is grid or not.

    Parameters
    ----------
    inputs_obj : CLIparser
        An object containing the command line input flags, specifically the number of threads.
    configuration : ConfigParser
        An object containing the configuration settings, including the computing mode (grid or local).
    """

    def __init__(self, inputs_obj: CLIparser, configuration: ConfigParser):
        self.assignment = False
        # NOTE: I expect that we need to change how 'see' the configuration of either local or grid mode, but this depends on other modifications.
        if configuration["COMPUTING"]["compmode"] == "grid":
            self.assignment = True

        threads = inputs_obj.flags.threads

        if not self.assignment:
            # Ensure calculated CPUs are at least 1 for non-grid mode
            self.highcpu = max(1, min(int(threads - 2), 12))
            self.midcpu = max(1, min(threads // 2, 6))
            self.lowcpu = 1
        else:  # Grid mode
            self.highcpu = 12
            self.midcpu = 6
            self.lowcpu = 2
