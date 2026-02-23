"""
Utilities for checking and updating ViroConstrictor from conda channels.

This module provides helpers to fetch package metadata from the Anaconda
API, perform an update using `conda` or `mamba`, and re-execute the updated
`viroconstrictor` binary.

"""

import configparser
import json
import os
import shutil
import subprocess
import sys
from typing import Any, NoReturn
from urllib import request

import packaging.version

from ViroConstrictor import __prog__, __version__
from ViroConstrictor.logging import log
from ViroConstrictor.userprofile import AskPrompts
from urllib.parse import urlparse

api_url = f"https://api.anaconda.org/release/bioconda/{__prog__.lower()}/latest"


def fetch_online_metadata() -> dict[str, Any] | None:
    """
    Fetch online metadata for the package from the Anaconda API.

    The function performs a blocking HTTP request to the Anaconda API URL
    configured in `api_url` and returns the parsed JSON response.

    Returns
    -------
    dict[str, Any] | None
        Parsed JSON metadata returned by the Anaconda API, or ``None`` if the
        request failed. A warning is logged when a connection or parsing
        error occurs.
    """
    try:
        parsed = urlparse(api_url)
        if parsed.scheme not in ("http", "https"):
            log.warning(f"Refusing to open URL with unsupported scheme: {parsed.scheme}")
            return None

        online_metadata = request.urlopen(api_url, timeout=60)
    except Exception as e:
        log.warning("Unable to connect to Anaconda API\n" f"{e}")
        return None
    return json.loads(online_metadata.read().decode("utf-8"))


def post_install(sysargs: list[str], online_version: packaging.version.Version) -> NoReturn:
    """
    Notify the user about a successful update and re-execute the updated binary.

    Parameters
    ----------
    sysargs : list[str]
        Command-line argument vector to pass to the new `viroconstrictor`
        process when re-executing the program.
    online_version : packaging.version.Version
        The version that was installed.

    Raises
    ------
    SystemExit
        If the `viroconstrictor` executable cannot be located after the
        update, the function logs an error and exits the process.

    Notes
    -----
    The function uses :func:`os.execv` to replace the running process so the
    updated code is loaded into memory.
    """
    log.info(f"ViroConstrictor updated to version [bold yellow]{online_version}[/bold yellow]")

    # Find the viroconstrictor executable path in the updated environment
    viroconstrictor_path = shutil.which("viroconstrictor")
    if viroconstrictor_path is None:
        log.error("Could not find viroconstrictor executable after update")
        sys.exit(1)

    # Use os.execv to replace the current process with the new viroconstrictor process.
    # Ensure that argv[0] matches the executable path by rebuilding the argument list.
    os.execv(viroconstrictor_path, [viroconstrictor_path] + sysargs[1:])


def _get_online_version() -> packaging.version.Version | None:
    """
    Retrieve the latest package version from the Anaconda API metadata.

    The function queries the online metadata via :func:`fetch_online_metadata`
    and parses the version string from the first distribution entry if
    available.

    Returns
    -------
    packaging.version.Version | None
        A parsed :class:`packaging.version.Version` instance for the latest
        available release, or ``None`` when metadata could not be fetched or
        no version could be determined.
    """

    online_metadata = fetch_online_metadata()
    if online_metadata is None:
        return None
    if latest_online_release := online_metadata.get("distributions", [])[0]:
        if version := latest_online_release.get("version"):
            return packaging.version.parse(version)
    return None


def _get_conda_prefix() -> str | None:
    """
    Get the active Conda/Mamba environment prefix from the environment.

    Returns
    -------
    str | None
        The value of the ``CONDA_PREFIX`` environment variable if set,
        otherwise ``None``.
    """

    conda_prefix = os.environ.get("CONDA_PREFIX")
    return conda_prefix or None


def _build_update_cmd(online_version: packaging.version.Version, conda_prefix: str) -> list[str]:
    """
    Construct the Conda/Mamba command used to perform the package update.

    Parameters
    ----------
    online_version : packaging.version.Version
        The target version to install.
    conda_prefix : str
        Path to the active Conda/Mamba environment (``CONDA_PREFIX``).

    Returns
    -------
    list[str]
        A command line argument list suitable for passing to
        :func:`subprocess.run` to perform the environment install.
    """

    conda_exe = "mamba" if shutil.which("mamba") else "conda"
    return [
        conda_exe,
        "install",
        "-y",
        "-p",
        conda_prefix,
        "-c",
        "bioconda",
        "-c",
        "conda-forge",
        "--force-reinstall",
        f"{__prog__.lower()}={online_version}",
    ]


def _run_update(online_version: packaging.version.Version) -> bool:
    """
    Execute the update workflow to install the specified version.

    This function first attempts to uninstall any pip-installed instance of
    the program, then runs a Conda/Mamba install into the active environment
    determined by ``CONDA_PREFIX``.

    Parameters
    ----------
    online_version : packaging.version.Version
        The version to install.

    Returns
    -------
    bool
        ``True`` when the update command returned a successful exit code,
        otherwise ``False``.

    Notes
    -----
    Any errors from the Conda/Mamba command are logged; the function does
    not raise on failure but returns ``False`` to allow the caller to handle
    post-failure behaviour.
    """

    subprocess.run([sys.executable, "-m", "pip", "uninstall", "-y", __prog__.lower()], capture_output=True, check=False)
    conda_prefix = _get_conda_prefix()
    if not conda_prefix:
        log.error(
            "ViroConstrictor update requires an active Conda or Mamba environment "
            "(environment variable CONDA_PREFIX is not set). "
            "Please activate the ViroConstrictor environment and retry the update."
        )
        return False

    update_cmd = _build_update_cmd(online_version, conda_prefix)
    result = subprocess.run(update_cmd, capture_output=True, text=True, check=False)

    if result.returncode != 0:
        if result.stderr:
            log.error("Conda/mamba reported the following error during the ViroConstrictor update:\n" f"{result.stderr.strip()}")
        else:
            log.error("Conda/mamba update command failed but did not produce any error output.")
    return result.returncode == 0


def _handle_update_result(sysargs: list[str], online_version: packaging.version.Version, update_successful: bool) -> None:
    """
    Handle the result of an attempted update.

    Parameters
    ----------
    sysargs : list[str]
        Command-line arguments passed to the current process; forwarded to
        :func:`post_install` when the update succeeded.
    online_version : packaging.version.Version
        The attempted target version.
    update_successful : bool
        ``True`` if the update command reported success.

    Returns
    -------
    None

    Notes
    -----
    On success the function delegates to :func:`post_install` which
    replaces the running process; on failure it logs an error describing
    possible transient solver issues.
    """

    if update_successful:
        post_install(sysargs, online_version)
        return
    log.error(
        f"ViroConstrictor update process to version [bold yellow]{online_version}[/bold yellow] failed at package solver stage.\n"
        "This might indicate that a new version of ViroConstrictor is uploaded to bioconda but download-servers are not ready yet.\n"
        "Please try again later."
    )


def _update_if_available(sysargs: list[str], local_version: packaging.version.Version) -> None:
    """
    Check for an available update and run it automatically if newer.

    Parameters
    ----------
    sysargs : list[str]
        Command-line arguments forwarded to the re-executed binary on success.
    local_version : packaging.version.Version
        Currently installed local version.

    Returns
    -------
    None

    Notes
    -----
    If an online version is found and is newer than ``local_version``, the
    function attempts to perform the update and will call
    :func:`post_install` on success. Errors are logged and the function
    returns without modifying the running process.
    """

    online_version = _get_online_version()
    if not online_version or local_version >= online_version:
        return

    log.info(f"Updating ViroConstrictor to latest version: [bold yellow]{online_version}[/bold yellow]")
    try:
        update_successful = _run_update(online_version)
    except Exception as exc:
        log.error(
            f"ViroConstrictor update process to version [bold yellow]{online_version}[/bold yellow] failed at package solver stage.\n"
            "This might indicate that a new version of ViroConstrictor is uploaded to bioconda but download-servers are not ready yet.\n"
            "Please try again later.\n"
            f"[red]Exception details:[/red] {exc}"
        )
        log.warning(f"Continuing with current version: [bold red]{local_version}[/bold red]")
        return

    if update_successful:
        post_install(sysargs, online_version)
        return

    _handle_update_result(sysargs, online_version, update_successful)
    log.warning(f"Continuing with current version: [bold red]{local_version}[/bold red]")


def _prompt_for_update(local_version: packaging.version.Version) -> packaging.version.Version | None:
    """
    Prompt the user to confirm updating to a newer online version.

    The function fetches the latest online version and, if it is newer
    than ``local_version``, displays a prompt using :class:`AskPrompts` to
    request consent. When the user consents the parsed online
    :class:`packaging.version.Version` is returned.

    Parameters
    ----------
    local_version : packaging.version.Version
        Currently installed package version.

    Returns
    -------
    packaging.version.Version | None
        The online version to update to when the user agrees, or ``None``
        when no update is available or the user declines.
    """

    online_version = _get_online_version()
    if not online_version or local_version >= online_version:
        return None

    if (
        AskPrompts(
            f"""
There's a new version of ViroConstrictor available.

Current version: [bold red]{local_version}[/bold red]
Latest version: [bold green]{online_version}[/bold green]\n""",
            "Do you want to update? [yes/no] ",
            ["yes", "no"],
            fixedchoices=True,
        )
        == "yes"
    ):
        return online_version

    log.info(f"Skipping update to version: [bold yellow]{online_version}[/bold yellow]")
    return None


def update(sysargs: list[str], conf: configparser.ConfigParser) -> None:
    """
    Check for and optionally install updates for ViroConstrictor.

    The function consults the local package version and the Anaconda API to
    determine whether an update is available. Behavior is controlled by the
    provided configuration parser in the `GENERAL` section:

    - If ``auto_update`` is ``yes``, the update is attempted automatically.
    - If ``auto_update`` is ``no`` and ``ask_for_update`` is ``yes``, the user
        is prompted before updating.

    Parameters
    ----------
    sysargs : list[str]
        Command line arguments passed to the current process; used when
        re-executing the program after a successful update.
    conf : configparser.ConfigParser
        Configuration with a `GENERAL` section that may contain
        ``auto_update`` and ``ask_for_update`` keys.

    Returns
    -------
    None

    Notes
    -----
    On success the function calls :func:`post_install` which replaces the
    running process with the updated ``viroconstrictor`` executable. Errors
    during the update process are logged and the function returns without
    changing the running process.
    """
    local_version = packaging.version.parse(__version__)

    autocontinue = conf.getboolean("GENERAL", "auto_update", fallback=False)
    ask_prompt = False if autocontinue else conf.getboolean("GENERAL", "ask_for_update", fallback=False)

    if autocontinue:
        _update_if_available(sysargs, local_version)
        return

    if not ask_prompt:
        return

    online_version = _prompt_for_update(local_version)
    if not online_version:
        return

    log.info(f"Updating ViroConstrictor to latest version: [bold yellow]{online_version}[/bold yellow]")
    try:
        update_successful = _run_update(online_version)
    except Exception:
        log.error(
            f"ViroConstrictor update process to version [bold yellow]{online_version}[/bold yellow] failed at package solver stage.\n"
            "This might indicate that a new version of ViroConstrictor is uploaded to bioconda but download-servers are not ready yet.\n"
            "Please try again later."
        )
        log.warning(f"Continuing with current version: [bold red]{local_version}[/bold red]")
        return

    _handle_update_result(sysargs, online_version, update_successful)
    if update_successful:
        return
    log.warning(f"Continuing with current version: [bold red]{local_version}[/bold red]")
