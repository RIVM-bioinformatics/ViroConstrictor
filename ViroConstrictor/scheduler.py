"""Module for the Scheduler enum."""

import os
import shutil
from configparser import ConfigParser
from enum import Enum
from logging import Logger


class Scheduler(Enum):
    """
    Enum to represent the scheduler type.
    The first entry in the tuples is the name of the executor plugin (except for auto),
    the rest are optional aliases.
    """

    LOCAL = ("local", "none", "")
    SLURM = ("slurm",)
    LSF = ("lsf",)
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
        return any(scheduler_str.lower() in scheduler.value for scheduler in cls)

    @classmethod
    def determine_scheduler(
        cls, scheduler_str: str, user_config: ConfigParser, log: Logger
    ) -> "Scheduler":
        """
        Determine the scheduler type from a string.
        It handles four options:
        1. If a string is provided, it converts it to a Scheduler enum member.
        2. It checks the user configuration for the scheduler. (if provided)
        3. It checks the environment variables for the scheduler.
        4. If no string or environment variables is found, it checks using drmaa (if available).
        If no scheduler is found, it returns Scheduler.LOCAL.
        """

        def _scheduler_by_env() -> Scheduler | None:
            """
            Determine the scheduler type from the environment variable.
            """
            if shutil.which("sbatch") or "SLURM_JOB_ID" in os.environ:
                return Scheduler.SLURM
            if shutil.which("bsub") or "LSB_JOBID" in os.environ:
                return Scheduler.LSF
            else:  # This should not return Scheduler.NONE, because then the 3rd option in determine_scheduler would never be used.
                return None

        log.debug("Determining scheduler...")
        scheduler: Scheduler | None

        # 1. From argument
        if scheduler_str:
            if not Scheduler.is_valid(scheduler_str):
                log.warning(
                    "Invalid scheduler string: '%s', using non-grid mode", scheduler_str
                )
                return Scheduler.LOCAL
            log.debug("Scheduler determined from string: '%s'", scheduler_str)
            scheduler = Scheduler.from_string(scheduler_str)
            if scheduler == Scheduler.AUTO:
                log.debug("Scheduler set to AUTO, trying to determine automatically")
            else:
                log.debug("Scheduler set to '%s'", scheduler.name)
                return scheduler

        # 2. From user configuration
        config_scheduler = user_config["COMPUTING"].get("scheduler", "")
        if config_scheduler:
            if not Scheduler.is_valid(config_scheduler):
                log.warning(
                    "Invalid scheduler in config: '%s', using non-grid mode",
                    config_scheduler,
                )
                return Scheduler.LOCAL
            log.debug("Scheduler determined from config: '%s'", config_scheduler)
            return Scheduler.from_string(config_scheduler)

        # 3. from environment variables
        scheduler = _scheduler_by_env()
        if scheduler is not None:
            log.debug("Scheduler determined from environment: '%s'", scheduler.name)
            return scheduler

        # 4. from drmaa (if available)
        try:
            # The DRMAA library is a wrapper, and fails to import if the C library is not installed (libdrmaa.so).
            import drmaa  # pylint: disable=import-outside-toplevel

            with drmaa.Session() as session:
                scheduler_name = session.drmsInfo
                log.debug("Scheduler determined from DRMAA: '%s'", scheduler_name)
                return Scheduler.from_string(scheduler_name)

        except Exception as e:  # pylint: disable=broad-except
            log.debug(f"DRMAA not available: {e}")

        log.info(
            "[yellow]No scheduler detected, running in non-grid mode. Please check your configuration or environment variables.[/yellow]"
        )
        return Scheduler.LOCAL
