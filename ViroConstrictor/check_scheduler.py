"""Module for the Scheduler enum."""

from enum import Enum


class Scheduler(Enum):
    """
    Enum to represent the scheduler type.
    """

    NONE = ("none", "")
    SLURM = ("slurm",)
    LSF = ("lsf",)

    @classmethod
    def supported_schedulers(cls) -> list[str]:
        """Get a list of all scheduler names."""
        return [scheduler.name for scheduler in cls]

    @classmethod
    def from_string(cls, scheduler_str: str) -> "Scheduler":
        """Convert a string to a Scheduler enum member."""
        for scheduler in cls:
            if scheduler_str.lower() in scheduler.value:
                return scheduler
        raise ValueError(
            f"Invalid scheduler: {scheduler_str}. Must be one of {cls.supported_schedulers()}."
        )

    @classmethod
    def is_valid(cls, scheduler_str: str) -> bool:
        """Check if the given string is a valid scheduler."""
        return any(scheduler_str.lower() in scheduler.value for scheduler in cls)
