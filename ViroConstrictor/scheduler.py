"""Module for the Scheduler enum."""

import os
import shutil
from configparser import ConfigParser
from enum import Enum
from typing import Optional

from ViroConstrictor.logging import log


class Scheduler(Enum):
    """
    Enum to represent the scheduler type.
    The first entry in the tuples is the name of the executor plugin (except for auto),
    the rest are optional aliases.
    """

    LOCAL = ("local", "none", "")
    SLURM = ("slurm",)
    LSF = ("lsf",)
    DRYRUN = ("dryrun",)
    AUTO = ("auto",)

    @classmethod
    def supported_schedulers(cls) -> list[str]:
        """Get a list of all scheduler names."""
        return [scheduler.name for scheduler in cls]

    @classmethod
    def from_string(cls, scheduler_str: str) -> "Scheduler":
        """Convert a string to a Scheduler enum member."""
        for scheduler in cls:
            if scheduler_str.lower().strip() in scheduler.value:
                return scheduler
        raise ValueError(
            f"Invalid scheduler: {scheduler_str}. Must be one of {cls.supported_schedulers()}."
        )

    @classmethod
    def is_valid(cls, scheduler_str: str) -> bool:
        """Check if the given string is a valid scheduler."""
        return any(
            scheduler_str.lower().strip() in scheduler.value for scheduler in cls
        )

    @classmethod
    def _scheduler_from_argument(
        cls,
        scheduler_str: str,
    ) -> Optional["Scheduler"]:
        if scheduler_str:
            if not cls.is_valid(scheduler_str):
                log.warning(
                    f"Invalid scheduler string: '{scheduler_str}', using non-grid mode"
                )
                return cls.LOCAL
            scheduler = cls.from_string(scheduler_str)
            if scheduler == cls.AUTO:
                log.debug(
                    "Helper functionality :: Scheduler :: Scheduler set to AUTO, trying to determine automatically"
                )
                return None
            log.debug(
                f"Helper functionality :: Scheduler :: Scheduler determined from string: '{scheduler_str}'"
            )
            return scheduler
        return None

    @classmethod
    def _scheduler_from_config(cls, user_config: ConfigParser) -> Optional["Scheduler"]:
        config_exec_mode = user_config["COMPUTING"].get("compmode", "")
        if config_exec_mode.lower() == "local":
            log.debug(
                "Helper functionality :: Scheduler :: Execution mode set to local in config, using LOCAL scheduler."
            )
            return cls.LOCAL

        config_scheduler = user_config["COMPUTING"].get("scheduler", "")
        if config_scheduler:
            if not cls.is_valid(config_scheduler):
                log.warning(
                    f"Invalid scheduler in config: '{config_scheduler}', using non-grid mode"
                )
                return cls.LOCAL
            log.debug(
                f"Helper functionality :: Scheduler :: Scheduler determined from config: '{config_scheduler}'"
            )
            return cls.from_string(config_scheduler)
        return None

    @classmethod
    def _scheduler_from_environment(cls) -> Optional["Scheduler"]:
        if shutil.which("sbatch") or "SLURM_JOB_ID" in os.environ:
            log.debug(
                "Helper functionality :: Scheduler :: Scheduler found in environment: SLURM"
            )
            return cls.SLURM
        if shutil.which("bsub") or "LSB_JOBID" in os.environ:
            log.debug(
                "Helper functionality :: Scheduler :: Scheduler found in environment: LSF"
            )
            return cls.LSF
        log.debug(
            "Helper functionality :: Scheduler :: No scheduler found in environment variables or executables."
        )
        return None

    @classmethod
    def _scheduler_from_drmaa(cls) -> Optional["Scheduler"]:
        try:
            import drmaa  # pylint: disable=import-outside-toplevel

            with drmaa.Session() as session:
                scheduler_name = session.drmsInfo
                log.debug(
                    f"Helper functionality :: Scheduler :: Scheduler determined from DRMAA: '{scheduler_name}'"
                )
                return cls.from_string(scheduler_name)
        except RuntimeError as e:
            log.debug(f"Helper functionality :: Scheduler :: DRMAA not available: {e}")
        return None

    @classmethod
    def determine_scheduler(
        cls,
        scheduler_str: str,
        user_config: ConfigParser,
        dryrun_arg: bool,
    ) -> "Scheduler":
        """Determine the scheduler type from argument, config, env, or DRMAA."""

        log.debug("Helper functionality :: Scheduler :: Determining scheduler...")
        if dryrun_arg:
            log.debug(
                "Helper functionality :: Scheduler :: Dry-run mode set on the commandline, using DRYRUN scheduler."
            )
            return cls.DRYRUN

        if scheduler_str:
            scheduler = cls._scheduler_from_argument(scheduler_str)
            if scheduler is not None:
                log.debug(
                    f"Helper functionality :: Scheduler :: Scheduler selected from argument: '{scheduler.name}'"
                )
                return scheduler

        if user_config:
            scheduler = cls._scheduler_from_config(user_config)
            if scheduler is not None:
                log.debug(
                    f"Helper functionality :: Scheduler :: Scheduler selected from config: '{scheduler.name}'"
                )
                return scheduler

        scheduler = cls._scheduler_from_environment()
        if scheduler is not None:
            log.debug(
                f"Helper functionality :: Scheduler :: Scheduler selected from environment: '{scheduler.name}'"
            )
            return scheduler

        scheduler = cls._scheduler_from_drmaa()
        if scheduler is not None:
            log.debug(
                f"Helper functionality :: Scheduler :: Scheduler selected from DRMAA: '{scheduler.name}'"
            )
            return scheduler

        log.info(
            "[yellow]No scheduler detected, running in non-grid mode. "
            "Please check your configuration or environment variables.[/yellow]"
        )
        return cls.LOCAL
