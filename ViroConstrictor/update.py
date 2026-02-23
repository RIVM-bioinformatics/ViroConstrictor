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

repo_channels = ("bioconda", "conda-forge")
api_url = f"https://api.anaconda.org/release/bioconda/{__prog__.lower()}/latest"

# TODO: current implementation contains a lot of code duplication. We should refactor this to avoid that but it is acceptable for the current release. Target release for fixing this is 1.7.0

# TODO: Current implementation is too complex with one huge updating method. This should be fixed in a future release, but acceptable for now. Target release for fixing this is 1.7.0


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
        online_metadata = request.urlopen(api_url)
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


# TODO: split this into smaller functions
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
    online_version = None

    autocontinue = conf["GENERAL"]["auto_update"] == "yes"

    ask_prompt = None
    if autocontinue == False:
        ask_prompt = conf["GENERAL"]["ask_for_update"] == "yes"

    if autocontinue:
        online_metadata = fetch_online_metadata()
        if online_metadata is None:
            return
        if latest_online_release := online_metadata.get("distributions", [])[0]:
            release_metadata = latest_online_release
            online_version = packaging.version.parse(release_metadata.get("version"))

            if online_version is not None and (local_version < online_version):
                log.info(f"Updating ViroConstrictor to latest version: [bold yellow]{online_version}[/bold yellow]")
                update_successful = False
                try:
                    # Always use mamba/conda for updates
                    # First, uninstall any pip-installed version to avoid conflicts
                    subprocess.run([sys.executable, "-m", "pip", "uninstall", "-y", __prog__.lower()], capture_output=True, check=False)

                    # Determine whether to use mamba or conda
                    conda_exe = "mamba" if shutil.which("mamba") else "conda"

                    # Determine target conda environment prefix
                    conda_prefix = os.environ.get("CONDA_PREFIX")
                    if not conda_prefix:
                        log.error(
                            "ViroConstrictor update requires an active Conda or Mamba environment "
                            "(environment variable CONDA_PREFIX is not set). "
                            "Please activate the ViroConstrictor environment and retry the update."
                        )
                        log.warning(f"Continuing with current version: [bold red]{local_version}[/bold red]")
                        return

                    # Build the install command (use install, not update, to force reinstallation)
                    update_cmd = [
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

                    result = subprocess.run(update_cmd, capture_output=True, text=True, check=False)
                    if result.returncode != 0:
                        # Surface solver/installation errors from conda/mamba to help debugging
                        if result.stderr:
                            log.error(
                                "Conda/mamba reported the following error during the ViroConstrictor update:\n"
                                f"{result.stderr.strip()}"
                            )
                        else:
                            log.error("Conda/mamba update command failed but did not produce any error output.")

                    update_successful = result.returncode == 0
                except Exception as e:
                    log.error(
                        f"ViroConstrictor update process to version [bold yellow]{online_version}[/bold yellow] failed at package solver stage.\n"
                        f"This might indicate that a new version of ViroConstrictor is uploaded to bioconda but download-servers are not ready yet.\n"
                        f"Please try again later.\n[red]Exception details:[/red] {e}"
                    )

                    log.warning(f"Continuing with current version: [bold red]{local_version}[/bold red]")
                    return
                if update_successful:
                    post_install(sysargs, online_version)
                else:
                    log.error(
                        f"ViroConstrictor update process to version [bold yellow]{online_version}[/bold yellow] failed at package solver stage.\nThis might indicate that a new version of ViroConstrictor is uploaded to bioconda but download-servers are not ready yet.\nPlease try again later."
                    )
                    log.warning(f"Continuing with current version: [bold red]{local_version}[/bold red]")
                return
    if not ask_prompt:
        return

    online_metadata = fetch_online_metadata()
    if online_metadata is None:
        return
    if latest_online_release := online_metadata.get("distributions", [])[0]:
        online_version = packaging.version.parse(latest_online_release.get("version"))

        if local_version < online_version:
            if (
                AskPrompts(
                    f"""
There's a new version of ViroConstrictor available.

Current version: [bold red]{local_version}[/bold red]
Latest version: [bold green]{online_version}[/bold green]\n""",
                    "Do you want to update? \[yes/no] ",
                    ["yes", "no"],
                    fixedchoices=True,
                )
                == "yes"
            ):
                log.info(f"Updating ViroConstrictor to latest version: [bold yellow]{online_version}[/bold yellow]")

                update_successful = False
                try:
                    # Always use mamba/conda for updates
                    # First, uninstall any pip-installed version to avoid conflicts
                    subprocess.run([sys.executable, "-m", "pip", "uninstall", "-y", __prog__.lower()], capture_output=True, check=False)

                    # Determine whether to use mamba or conda
                    conda_exe = "mamba" if shutil.which("mamba") else "conda"

                    # Determine target conda environment prefix
                    conda_prefix = os.environ.get("CONDA_PREFIX")
                    if not conda_prefix:
                        log.error(
                            "Unable to perform self-update: ViroConstrictor must be running inside a conda or mamba "
                            "environment so that CONDA_PREFIX is defined. Please activate the appropriate environment "
                            "and rerun ViroConstrictor with the update option."
                        )
                        log.warning(f"Skipping update; continuing with current version: [bold red]{local_version}[/bold red]")
                        return

                    # Build the install command (use install, not update, to force reinstallation)
                    update_cmd = [
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

                    result = subprocess.run(update_cmd, capture_output=True, text=True, check=False)

                    update_successful = result.returncode == 0
                except Exception:
                    update_successful = False

                    log.error(
                        f"ViroConstrictor update process to version [bold yellow]{online_version}[/bold yellow] failed at package solver stage.\nThis might indicate that a new version of ViroConstrictor is uploaded to bioconda but download-servers are not ready yet.\nPlease try again later."
                    )
                    log.warning(f"Continuing with current version: [bold red]{local_version}[/bold red]")
                if update_successful:
                    post_install(sysargs, online_version)
                else:
                    log.error(
                        f"ViroConstrictor update process to version [bold yellow]{online_version}[/bold yellow] failed at package solver stage.\nThis might indicate that a new version of ViroConstrictor is uploaded to bioconda but download-servers are not ready yet.\nPlease try again later."
                    )
                    log.warning(f"Continuing with current version: [bold red]{local_version}[/bold red]")
                return
            log.info(f"Skipping update to version: [bold yellow]{online_version}[/bold yellow]")
            return
        return
