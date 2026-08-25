"""Provide lightweight logging helpers for the container_manager package."""

import logging
import os

_STATE = {"configured": False}


def _configure_once() -> None:
    """Configure root logging once for this package."""
    if _STATE["configured"]:
        return

    level_name = os.environ.get("CONTAINER_MANAGER_LOG_LEVEL", "INFO").upper()
    level = getattr(logging, level_name, logging.INFO)
    logging.basicConfig(level=level, format="%(asctime)s %(levelname)s %(name)s: %(message)s")
    _STATE["configured"] = True


def get_logger(name: str) -> logging.Logger:
    """Return a logger configured for container_manager use."""
    _configure_once()
    return logging.getLogger(name)
