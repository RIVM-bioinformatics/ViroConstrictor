import configparser
import contextlib
import json
import os
import shutil
import subprocess
import sys
from typing import Any, NoReturn
from urllib import request

import packaging.version
from rich import print

from ViroConstrictor import __prog__, __version__
from ViroConstrictor.logging import log
from ViroConstrictor.userprofile import AskPrompts

repo_channels = ("bioconda", "conda-forge")
api_url = f"https://api.anaconda.org/release/bioconda/{__prog__.lower()}/latest"


@contextlib.contextmanager
def silence_stdout_stderr():
    """
    This function temporarily redirects stdout and stderr to a null device to silence them.
    """
    sderr_fd = sys.stderr.fileno()
    sdout_fd = sys.stdout.fileno()
    orig_stderr_fd = os.dup(sderr_fd)
    orig_stdout_fd = os.dup(sdout_fd)
    null_fd = os.open(os.devnull, os.O_WRONLY)
    os.dup2(null_fd, sdout_fd)
    os.dup2(null_fd, sderr_fd)
    try:
        yield
    finally:
        os.dup2(orig_stdout_fd, sdout_fd)
        os.dup2(orig_stderr_fd, sderr_fd)
        os.close(null_fd)
        os.close(orig_stderr_fd)
        os.close(orig_stdout_fd)


def fetch_online_metadata() -> dict[str, Any] | None:
    """This function fetches online metadata from the Anaconda API URL and returns it as a dictionary, or returns
    None if there is an exception.

    Returns
    -------
        A dictionary with string keys and values of any type, or `None` if an exception occurs while trying to connect to the Anaconda API.

    """
    try:
        online_metadata = request.urlopen(api_url)
    except Exception as e:
        log.warning("Unable to connect to Anaconda API\n" f"{e}")
        return None
    return json.loads(online_metadata.read().decode("utf-8"))


def post_install(sysargs: list[str], online_version: packaging.version.Version) -> NoReturn:
    """This function prints a message indicating the updated version of ViroConstrictor and replaces
    the current process with a new one running the updated ViroConstrictor command.

    Parameters
    ----------
    sysargs : list[str]
        A list of command-line arguments to be passed to the new process.
    online_version : LooseVersion
        The version number of the updated ViroConstrictor package that is available online.

    """
    log.info(f"ViroConstrictor updated to version [bold yellow]{online_version}[/bold yellow]")
    
    # Find the viroconstrictor executable path in the updated environment
    viroconstrictor_path = shutil.which("viroconstrictor")
    if viroconstrictor_path is None:
        log.error("Could not find viroconstrictor executable after update")
        sys.exit(1)
    
    # Use os.execv to replace the current process with the new viroconstrictor process
    # This ensures the updated code is loaded
    os.execv(viroconstrictor_path, sysargs)


# TODO: split this into smaller functions
def update(sysargs: list[str], conf: configparser.ConfigParser) -> None:
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
                update_succesfull = False
                try:
                    # Always use mamba/conda for updates
                    # First, uninstall any pip-installed version to avoid conflicts
                    subprocess.run(
                        [sys.executable, "-m", "pip", "uninstall", "-y", __prog__.lower()],
                        capture_output=True,
                        check=False
                    )
                    
                    # Determine whether to use mamba or conda
                    conda_exe = "mamba" if shutil.which("mamba") else "conda"
                    
                    # Build the install command (use install, not update, to force reinstallation)
                    update_cmd = [
                        conda_exe,
                        "install",
                        "-y",
                        "-p", os.environ["CONDA_PREFIX"],
                        "-c", "bioconda",
                        "-c", "conda-forge",
                        "--force-reinstall",
                        f"{__prog__.lower()}={online_version}"
                    ]
                    
                    result = subprocess.run(
                        update_cmd,
                        capture_output=True,
                        text=True,
                        check=False
                    )
                    
                    update_succesfull = result.returncode == 0
                except Exception as e:
                    log.error(
                        f"ViroConstrictor update process to version [bold yellow]{online_version}[/bold yellow] failed at package solver stage.\nThis might indicate that a new version of ViroConstrictor is uploaded to bioconda but download-servers are not ready yet.\nPlease try again later."
                    )

                    log.warning(f"Continuing with current version: [bold red]{local_version}[/bold red]")
                    return
                if update_succesfull:
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

                update_succesfull = False
                try:
                    # Always use mamba/conda for updates
                    # First, uninstall any pip-installed version to avoid conflicts
                    subprocess.run(
                        [sys.executable, "-m", "pip", "uninstall", "-y", __prog__.lower()],
                        capture_output=True,
                        check=False
                    )
                    
                    # Determine whether to use mamba or conda
                    conda_exe = "mamba" if shutil.which("mamba") else "conda"
                    
                    # Build the install command (use install, not update, to force reinstallation)
                    update_cmd = [
                        conda_exe,
                        "install",
                        "-y",
                        "-p", os.environ["CONDA_PREFIX"],
                        "-c", "bioconda",
                        "-c", "conda-forge",
                        "--force-reinstall",
                        f"{__prog__.lower()}={online_version}"
                    ]
                    
                    result = subprocess.run(
                        update_cmd,
                        capture_output=True,
                        text=True,
                        check=False
                    )
                    
                    update_succesfull = result.returncode == 0
                except Exception as e:
                    update_succesfull = False

                    log.error(
                        f"ViroConstrictor update process to version [bold yellow]{online_version}[/bold yellow] failed at package solver stage.\nThis might indicate that a new version of ViroConstrictor is uploaded to bioconda but download-servers are not ready yet.\nPlease try again later."
                    )
                    log.warning(f"Continuing with current version: [bold red]{local_version}[/bold red]")
                if update_succesfull:
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
