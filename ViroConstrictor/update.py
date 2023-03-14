import contextlib
import json
import os
import subprocess
import sys
from distutils.version import LooseVersion
from urllib import request

from mamba.api import install as mamba_install
from rich import print

from ViroConstrictor import __prog__, __version__

from .functions import color
from .userprofile import AskPrompts

repo_channels = ("bioconda", "conda-forge", "intel", "anaconda")


@contextlib.contextmanager
def silence_stdout_stderr():
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


def update(sysargs, conf):

    local_version = LooseVersion(__version__)
    online_version = None

    autocontinue = conf["GENERAL"]["auto_update"] == "yes"
    ask_prompt = not autocontinue and conf["GENERAL"]["ask_for_update"] == "yes"

    if autocontinue:
        try:
            online_metadata = request.urlopen(
                f"https://api.anaconda.org/release/bioconda/{__prog__.lower()}/latest"
            )
        except Exception as e:
            sys.stderr.write("Unable to connect to Anaconda API\n" f"{e}")
            return

        online_metadata = json.loads(online_metadata.read().decode("utf-8"))
        if latest_online_release := online_metadata.get("distributions")[0]:
            release_metadata = latest_online_release
            online_version = LooseVersion(release_metadata.get("version"))

            if online_version is not None and (local_version < online_version):
                print(
                    f"Updating ViroConstrictor to latest version: [bold yellow]{online_version}[/bold yellow]"
                )
                result = False
                with silence_stdout_stderr():
                    result = mamba_install(
                        os.environ["CONDA_PREFIX"],
                        (f"{__prog__.lower()} {online_version}",),
                        repo_channels,
                    )
                if result:
                    print(
                        f"ViroConstrictor updated to version [bold yellow]{online_version}[/bold yellow]]"
                    )
                    subprocess.run(sysargs)
                    sys.exit(0)
                print(
                    f"Failed to update ViroConstrictor to version [bold yellow]{online_version}[/bold yellow]]]"
                )
                print("Please update manually")
                print(
                    f"Continuing with current version: [bold red]{local_version}[/bold red]"
                )
                return
    if not ask_prompt:
        return

    try:
        online_metadata = request.urlopen(
            f"https://api.anaconda.org/release/bioconda/{__prog__.lower()}/latest"
        )
    except Exception as e:
        sys.stderr.write("Unable to connect to Anaconda API\n" f"{e}")
        return

    online_metadata = json.loads(online_metadata.read().decode("utf-8"))
    if latest_online_release := online_metadata.get("distributions")[0]:
        online_version = LooseVersion(latest_online_release.get("version"))

        if local_version < online_version:
            if (
                AskPrompts(
                    f"""
There's a new version of ViroConstrictor available.

Current version: {color.BOLD + color.RED}{local_version}{color.END}
Latest version: {color.BOLD + color.GREEN}{online_version}{color.END}\n""",
                    "Do you want to update? [yes/no]",
                    ["yes", "no"],
                    fixedchoices=True,
                )
                == "yes"
            ):

                print(
                    f"Updating ViroConstrictor to latest version: [bold yellow]{online_version}[/bold yellow]"
                )

                result = False
                with silence_stdout_stderr():
                    result = mamba_install(
                        os.environ["CONDA_PREFIX"],
                        (f"{__prog__.lower()} {online_version}",),
                        repo_channels,
                    )
                if result:
                    print(
                        f"ViroConstrictor updated to version [bold yellow]{online_version}[/bold yellow]"
                    )
                    subprocess.run(sysargs)
                    sys.exit(0)
                print(
                    f"Failed to update ViroConstrictor to version [bold yellow]{online_version}[/bold yellow]"
                )
                print("Please update manually")
                print(
                    f"Continuing with current version: [bold red]{local_version}[/bold red]"
                )
                return
            print(
                f"Skipping update to version: [bold yellow]{online_version}[/bold yellow]"
            )
            print("Continuing...")
            return
        return
