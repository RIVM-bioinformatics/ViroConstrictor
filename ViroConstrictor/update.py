import json
import subprocess
import sys
from distutils.version import LooseVersion
from urllib import request

from .functions import color
from .userprofile import AskPrompts
from .version import __version__


def update(args, conf):

    autocontinue = False
    ask_prompt = False

    if conf["GENERAL"]["auto_update"] == "yes":
        autocontinue = True

    if autocontinue is False:
        if conf["GENERAL"]["ask_for_update"] == "yes":
            ask_prompt = True

    if autocontinue is True:
        try:
            latest_release = request.urlopen(
                "https://api.github.com/repos/RIVM-bioinformatics/Viroconstrictor/releases"
            )
        except Exception as e:
            sys.stderr.write("Unable to connect to GitHub API\n" f"{e}")
            return

        latest_release = json.loads(latest_release.read().decode("utf-8"))[0]

        latest_release_tag = latest_release["tag_name"]
        latest_release_tag_tidied = LooseVersion(
            latest_release["tag_name"].lstrip("v").strip()
        )

        localversion = LooseVersion(__version__)

        if localversion < latest_release_tag_tidied:
            subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "pip",
                    "install",
                    "--upgrade",
                    f"git+https://github.com/RIVM-bioinformatics/ViroConstrictor@{latest_release_tag}",
                ],
                check=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )

            print(
                f"ViroConstrictor updated to {color.YELLOW + color.BOLD}{latest_release_tag}{color.END}"
            )

            subprocess.run(args)
            sys.exit(0)
        return

    if autocontinue is False and ask_prompt is False:
        return

    if autocontinue is False and ask_prompt is True:
        try:
            latest_release = request.urlopen(
                "https://api.github.com/repos/RIVM-bioinformatics/Viroconstrictor/releases"
            )
        except Exception as e:
            sys.stderr.write("Unable to connect to GitHub API\n" f"{e}")
            return

        latest_release = json.loads(latest_release.read().decode("utf-8"))[0]

        latest_release_tag = latest_release["tag_name"]
        latest_release_tag_tidied = LooseVersion(
            latest_release["tag_name"].lstrip("v").strip()
        )

        localversion = LooseVersion(__version__)

        if localversion < latest_release_tag_tidied:
            if (
                AskPrompts(
                    f"""
There's a new version of ViroConstrictor available.

Current version: {color.RED + color.BOLD}{'v' + __version__}{color.END}
Latest version: {color.GREEN + color.BOLD}{latest_release_tag}{color.END}\n""",
                    f"""Do you want to update? [yes/no] """,
                    ["yes", "no"],
                    fixedchoices=True,
                )
                == "yes"
            ):
                subprocess.run(
                    [
                        sys.executable,
                        "-m",
                        "pip",
                        "install",
                        "--upgrade",
                        f"git+https://github.com/RIVM-bioinformatics/ViroConstrictor@{latest_release_tag}",
                    ],
                    check=True,
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.DEVNULL,
                )

                print(
                    f"ViroConstrictor updated to {color.YELLOW + color.BOLD}{latest_release_tag}{color.END}"
                )

                subprocess.run(args)
                sys.exit(0)
            print(f"Skipping update to version {latest_release_tag}")
            print(f"Continuing...")
            return
        return
