import json
import os
import subprocess
import sys
from distutils.version import LooseVersion
from urllib import request

from ViroConstrictor import __version__

from .functions import color
from .userprofile import AskPrompts


def update(sysargs, conf):

    autocontinue = conf["GENERAL"]["auto_update"] == "yes"
    ask_prompt = not autocontinue and conf["GENERAL"]["ask_for_update"] == "yes"
    if autocontinue:
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

        if (
            localversion < latest_release_tag_tidied
            and localversion.version[0] == latest_release_tag_tidied.version[0]
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

            subprocess.run(sysargs)
            sys.exit(0)
        if (
            localversion < latest_release_tag_tidied
            and localversion.version[0] != latest_release_tag_tidied.version[0]
        ):
            if (
                AskPrompts(
                    f"""{color.RED + color.BOLD}There's a new version of ViroConstrictor available. This new version is a {color.UNDERLINE}major{color.END + color.RED + color.BOLD} update and cannot be installed automatically.{color.END}

Current version: {color.RED + color.BOLD}{'v' + __version__}{color.END}
Latest version: {color.GREEN + color.BOLD}{latest_release_tag}{color.END}

The auto-updater can't install major version changes for you, as this would (probably) severely break your installation.
If you wish to update to the newest version you will have to do so manually.

If you want to run ViroConstrictor with the current version then please turn off the auto-updater. 
If you won't turn off the auto-updater we'll keep nagging you about this until you manually updated to the newest version.
""",
                    "Do you want to turn off the auto-updater so you wont get this message again? (yes/no) ",
                    ("yes", "no"),
                    fixedchoices=True,
                )
                == "yes"
            ):
                conf["GENERAL"]["auto_update"] = "no"
                conf["GENERAL"]["ask_for_update"] = "yes"

                with open(
                    os.path.expanduser("~/.ViroConstrictor_defaultprofile.ini"), "w"
                ) as f:
                    conf.write(f)
                print("The ViroConstrictor auto-updater is now turned off")
                print(
                    f"Please re-run the ViroConstrictor command to execute the workflow with the current version ({'v' + __version__}) or update manually to the newest version"
                )
                sys.exit(0)
            print(
                "ViroConstrictor is unable to update itself to the newest version as this is a major version change that cannot be installed automatically.\nPlease update manually and try again or turn-off the auto-updater "
            )
            print("Exiting...")
            sys.exit(1)
        return

    if not ask_prompt:
        return

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

    if (
            localversion < latest_release_tag_tidied
            and localversion.version[0] == latest_release_tag_tidied.version[0]
        ):
        if AskPrompts(f"""
There's a new version of ViroConstrictor available.

Current version: {color.RED + color.BOLD}{'v' + __version__}{color.END}
Latest version: {color.GREEN + color.BOLD}{latest_release_tag}{color.END}\n""", """Do you want to update? [yes/no] """, ["yes", "no"], fixedchoices=True) == "yes":
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

            subprocess.run(sysargs)
            sys.exit(0)
        print(f"Skipping update to version {latest_release_tag}")
        print("Continuing...")
        return
    if localversion < latest_release_tag_tidied:
        print(
            f"{color.RED}There's a new version of ViroConstrictor available. This new version is a {color.UNDERLINE}major{color.END + color.RED} update and cannot be installed automatically.{color.END}"
        )
        print(f"""Current version: {color.RED + color.BOLD}v{__version__}{color.END}\nLatest version: {color.GREEN + color.BOLD}{latest_release_tag}{color.END}""")

        print("Continuing without updating...\n")
    return
