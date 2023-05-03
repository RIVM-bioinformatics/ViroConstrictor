import configparser
import contextlib
import json
import os
import subprocess
import sys
from distutils.version import LooseVersion
from typing import Any, NoReturn
from urllib import request

from mamba.api import install as mamba_install
from rich import print

from ViroConstrictor import __prog__, __version__
from ViroConstrictor.logging import log
from ViroConstrictor.userprofile import AskPrompts

repo_channels = ("bioconda", "conda-forge", "intel", "anaconda")
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


def post_install(sysargs: list[str], online_version: LooseVersion) -> NoReturn:
    """This function prints a message indicating the updated version of ViroConstrictor and runs a
    subprocess of the original ViroConstrictor command before exiting the system.

    Parameters
    ----------
    sysargs : list[str]
        A list of command-line arguments to be passed to a subprocess that will be run after the print
    statement.
    online_version : LooseVersion
        The version number of the updated ViroConstrictor package that is available online.

    """
    print(
        f"ViroConstrictor updated to version [bold yellow]{online_version}[/bold yellow]"
    )
    subprocess.run(sysargs)
    sys.exit(0)


# TODO: split this into smaller functions
def update(sysargs: list[str], conf: configparser.ConfigParser) -> None:
    local_version = LooseVersion(__version__)
    online_version = None

    autocontinue = conf["GENERAL"]["auto_update"] == "yes"
    ask_prompt = conf["GENERAL"]["ask_for_update"] == "yes"

    if autocontinue:
        online_metadata = fetch_online_metadata()
        if online_metadata is None:
            return
        if latest_online_release := online_metadata.get("distributions", [])[0]:
            release_metadata = latest_online_release
            online_version = LooseVersion(release_metadata.get("version"))

            if online_version is not None and (local_version < online_version):
                log.info(
                    f"Updating ViroConstrictor to latest version: [bold yellow]{online_version}[/bold yellow]"
                )
                update_succesfull = False
                try:
                    with silence_stdout_stderr():
                        update_succesfull = mamba_install(
                            os.environ["CONDA_PREFIX"],
                            (f"{__prog__.lower()} {online_version}",),
                            repo_channels,
                        )
                except Exception as e:
                    update_succesfull = False
                    #
                    log.error(
                        f"ViroConstrictpor update process to version [bold yellow]{online_version}[/bold yellow] failed at package solver stage.\nThis might indicate that a new version of ViroConstrictor is uploaded to bioconda but download-servers are not ready yet.\nPlease try again later."
                    )

                    log.warning(
                        f"Continuing with current version: [bold red]{local_version}[/bold red]"
                    )
                    return
                if update_succesfull:
                    post_install(sysargs, online_version)
                return
    if not ask_prompt:
        return

    online_metadata = fetch_online_metadata()
    if online_metadata is None:
        return
    if latest_online_release := online_metadata.get("distributions", [])[0]:
        online_version = LooseVersion(latest_online_release.get("version"))

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
                log.info(
                    f"Updating ViroConstrictor to latest version: [bold yellow]{online_version}[/bold yellow]"
                )

                update_succesfull = False
                try:
                    with silence_stdout_stderr():
                        update_succesfull = mamba_install(
                            os.environ["CONDA_PREFIX"],
                            (f"{__prog__.lower()} {online_version}",),
                            repo_channels,
                        )
                except Exception as e:
                    update_succesfull = False

                    log.error(
                        f"ViroConstrictpor update process to version [bold yellow]{online_version}[/bold yellow] failed at package solver stage.\nThis might indicate that a new version of ViroConstrictor is uploaded to bioconda but download-servers are not ready yet.\nPlease try again later."
                    )
                    log.warning(
                        f"Continuing with current version: [bold red]{local_version}[/bold red]"
                    )
                if update_succesfull:
                    post_install(sysargs, online_version)
                return
            log.info(
                f"Skipping update to version: [bold yellow]{online_version}[/bold yellow]"
            )
            return
        return
