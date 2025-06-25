# pylint: disable=C0103

"""
Construct and write configuration files for SnakeMake
"""

import multiprocessing
import os
from typing import Any

import yaml

from ViroConstrictor.check_scheduler import Scheduler
from ViroConstrictor.parser import CLIparser


def WriteYaml(data: dict, filepath: str) -> str:
    """WriteYaml takes a dictionary and a filepath, and writes the given dictionary to the filepath as a yaml
    file

    Parameters
    ----------
    data : dict
        The data to be written to the file.
    filepath : str
        The path to the file you want to write to.

    Returns
    -------
    filepath : str
        The filepath

    """
    if not os.path.exists(os.path.dirname(filepath)):
        os.makedirs(os.path.dirname(filepath))
    with open(filepath, "w") as file:
        yaml.dump(data, file, default_flow_style=False)
    return filepath


class GetSnakemakeRunDetails:
    def __init__(
        self,
        inputs_obj: CLIparser,
        samplesheetfilename: str,
        scheduler: Scheduler,
        outdirOverride: str = "",
    ) -> None:
        self.inputs = inputs_obj
        self.samplesheetfilename = samplesheetfilename
        self.scheduler = scheduler
        self.outdirOverride = outdirOverride
        self._snakemake_run_config()
        self._snakemake_run_params()

    def _snakemake_run_config(self) -> dict[str, Any]:
        self.snakemake_run_conf = {}
        cores = self._set_cores(self.inputs.flags.threads)
        configuration = self.inputs.user_config
        compmode = configuration["COMPUTING"]["compmode"]
        reproduction_mode = configuration["REPRODUCTION"]["repro_method"]
        reproduction_cache_path = configuration["REPRODUCTION"]["container_cache_path"]

        base_config = {
            "latency-wait": 60,
            "dryrun": self.inputs.flags.dryrun,
            "jobname": "ViroConstrictor_{name}.jobid{jobid}",
            "use-conda": reproduction_mode == "conda",
            "use-singularity": reproduction_mode == "containers",
            "container_cache": reproduction_cache_path,
            "restart-times": 3,
            "keep-going": True,
            "printshellcmds": False,
            "scheduler": "greedy",
            "force-incomplete": True,
        }

        if compmode == "grid":
            queuename = configuration["COMPUTING"]["queuename"]
            drmaa_cmds = {
                Scheduler.LSF: f' -q {queuename} -n {{threads}} -R "span[hosts=1]" -M {{resources.mem_mb}}',
                Scheduler.SLURM: f" -p {queuename} -c {{threads}} --mem={{resources.mem_mb}}",
            }
            base_config.update(
                {
                    "cores": 300,
                    "drmaa": drmaa_cmds[self.scheduler],
                    "drmaa-log-dir": "logs/drmaa",
                }
            )
            self.snakemake_run_conf = base_config
            return self.snakemake_run_conf
        base_config.update(
            {
                "cores": cores,
            }
        )
        self.snakemake_run_conf = base_config
        return self.snakemake_run_conf

    def _snakemake_run_params(self) -> dict[str, Any]:
        self.snakemake_run_parameters = {}
        configuration = self.inputs.user_config
        threads_highcpu: int = min(int(self.inputs.flags.threads - 2), 12)
        threads_midcpu: int = min(self.inputs.flags.threads // 2, 6)
        threads_lowcpu: int = 1
        if configuration["COMPUTING"]["compmode"] == "grid":
            threads_highcpu = 12
            threads_midcpu = 6
            threads_lowcpu = 2
        self.snakemake_run_parameters = {
            "sample_sheet": WriteYaml(
                self.inputs.samples_dict,
                f"{self.inputs.workdir}/{self.samplesheetfilename}.yaml",
            ),
            "computing_execution": configuration["COMPUTING"]["compmode"],
            "max_local_mem": self._get_max_local_mem(),
            "platform": self.inputs.flags.platform,
            "amplicon_type": self.inputs.flags.amplicon_type,
            "outdirOverride": self.outdirOverride,
            "threads": {
                "Alignments": threads_highcpu,
                "QC": threads_midcpu,
                "AdapterRemoval": threads_lowcpu,
                "PrimerRemoval": threads_highcpu,
                "Consensus": threads_midcpu,
                "Index": threads_lowcpu,
                "Typing": threads_lowcpu,
            },
        }
        return self.snakemake_run_parameters

    def _set_cores(self, cores: int) -> int:
        """Balance the requested number of cores to the available number of cores on the system to ensure analysis will not overload the system in local-analysis mode.
        If the number of cores requested is equal to the number of cores available, return the number of cores requested minus 2.
        If the number of cores requested is greater than the number of cores available, then return the
        number of cores available minus 2.
        Otherwise, return the number of cores requested

        Parameters
        ----------
        cores : int
            The number of cores to use for parallel processing given by the user.

        Returns
        -------
        cores : int
            The number of cores that can be used by ViroConstrictor during analysis.

        """
        available: int = multiprocessing.cpu_count()
        if cores == available:
            return cores - 2
        return available - 2 if cores > available else cores

    def _get_max_local_mem(self) -> int:
        """It returns the maximum amount of local memory that can be allocated to a single process on the
        current machine

        Returns
        -------
            The amount of memory available in MB.

        """
        avl_mem_bytes = os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES")
        return int(round(avl_mem_bytes / (1024.0**2) - 2000, -3))
