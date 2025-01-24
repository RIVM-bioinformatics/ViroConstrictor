# pylint: disable=C0103

"""
Read or write the configuration file for ViroConstrictor.
This is so a user doesn't have to give this information manually on every run
"""
import configparser
import os
import pathlib
import readline
import subprocess
import sys
from typing import Any

from rich import print
from rich.console import Console

from ViroConstrictor.functions import tabCompleter
from ViroConstrictor.logging import log
from ViroConstrictor.workflow.containers import containerization_installed


def FileExists(file: pathlib.Path) -> bool:
    """Return True if the file exists, False if it doesn't.

    Args:
        file: The file to check.

    Returns:
        True if the file exists, False if it doesn't.
    """
    return bool(os.path.isfile(file))


def FileIsPopulated(file: pathlib.Path) -> bool:
    """Checks if the given file is populated.

    Args:
        file: The file to check.

    Returns:
        True if the given file is populated, False otherwise.
    """
    return os.stat(file).st_size >= 1


def AskPrompts(
    intro: str,
    prompt: str,
    options: list,
    fixedchoices: bool = False,
    default: Any = None,
) -> str | None:
    """This function is used to ask the user a question and provide a list of options to choose from.
    A free-text user reply is also possible.

    The function takes 4 arguments:

    1. intro: This is the introduction text that will be displayed to the user.
    2. prompt: This is the question that will be asked to the user.
    3. options: This is a list of options that the user can choose from.
    4. fixedchoices: This is a boolean value that determines whether the user can enter a custom answer
    or not

    Parameters
    ----------
    intro
        This is the text that will be displayed before the prompt.
    prompt
        The prompt that will be displayed to the user
    options
        a list of options to choose from
    fixedchoices, optional
        If set to True, the user will be able to use the tab key to autocomplete the available options.

    Returns
    -------
        the reply variable.

    """
    completer = tabCompleter()
    if fixedchoices:
        completer.createListCompleter(options)
        readline.set_completer_delims("\t")
        readline.parse_and_bind("tab: complete")
        readline.set_completer(completer.listCompleter)
    else:
        readline.parse_and_bind("tab: complete")
        readline.set_completer(completer.pathCompleter)

    subprocess.call("/bin/clear", shell=False)

    print(f"[bold blue]{'='*75}[/bold blue]")
    print(intro)
    print(f"[bold blue]{'='*75}[/bold blue]")
    while "the answer is invalid":
        if fixedchoices:
            reply = Console(soft_wrap=True).input(prompt)
            if reply in options:
                return reply
            if reply == "quit":
                print("Quitting...")
                sys.exit(-1)
            else:
                print(
                    "The given answer was invalid. Please choose one of the available options\n"
                )
        if not fixedchoices:
            reply = input(prompt).strip()
            if reply == "quit":
                sys.exit(-1)
            # if reply is empty and a default value is given, return the default value
            return default if not reply and default else reply


def BuildConfig(file: pathlib.Path) -> None:
    """Function asks the user a series of questions and writes the answers to a config file

    Parameters
    ----------
    file
        The file to write the config to.

    """
    # pylint: disable=C0301
    if os.path.exists(file):
        os.remove(file)

    conf_object = configparser.ConfigParser()

    conf_object["COMPUTING"] = {  # type: ignore
        "compmode": AskPrompts(
            """ViroConstrictor can run in two computing-modes.
[yellow underline]local[/yellow underline] or [yellow underline]HPC/Grid[/yellow underline]
Please specify the computing-mode that you wish to use for ViroConstrictor.""",
            """Do you wish to run ViroConstrictor in [yellow]local[/yellow] or [yellow]grid[/yellow] mode? \[local/grid] """,
            ["local", "grid"],
            fixedchoices=True,
        )
    }

    if conf_object["COMPUTING"]["compmode"] == "grid":
        conf_object["COMPUTING"]["queuename"] = AskPrompts(  # type: ignore
            """Grid mode has been chosen. Please enter the name of computing-queue that you wish to use on your grid/HPC cluster.\nThis is necessary so ViroConstrictor will send all the various tasks to the correct (remote) computers.\n\n[bold underline yellow]Please note that this is case-sensitive[/bold underline yellow]\n""",
            "Please specify the name of the Queue on your grid/HPC cluster that you wish to use. [free text] ",
            [],
            fixedchoices=False,
        )

    conf_object["GENERAL"] = {  # type: ignore
        "auto_update": AskPrompts(
            """ViroConstrictor can check and update itself everytime you run it.
Please specify whether you wish to enable the auto-update feature.""",
            """Do you wish to enable the auto-update feature? \[yes/no] """,
            ["yes", "no"],
            fixedchoices=True,
        )
    }

    if conf_object["GENERAL"]["auto_update"] == "no":
        conf_object["GENERAL"]["ask_for_update"] = AskPrompts(  # type: ignore
            """ViroConstrictor will not automatically update itself, but ViroConstrictor can still check for updates and ask you if you wish to update.""",
            """Do you want ViroConstrictor to [yellow underline]ask you[/yellow underline] to update everytime a new update is available? \[yes/no] """,
            ["yes", "no"],
            fixedchoices=True,
        )

    if containerization_installed is False:
        conf_object["REPRODUCTION"] = {
            "repro_method": "conda",
            "container_cache_path": None,
        }
    else:
        conf_object["REPRODUCTION"] = {"repro_method": "containers"}
        conf_object["REPRODUCTION"]["container_cache_path"] = AskPrompts(
            f"""ViroConstrictor will use containers to run the various analysis steps in a reproducible manner.\nHowever, to speed up the workflow and to allow offline-use, ViroConstrictor will cache the containers on your local machine.\nThis directory will be used to locally store the containers.\nIf you do not provide a path, the default path will be used. ({str(pathlib.Path.home())}/.viroconstrictor/containers)""",
            """Please specify the path to the container cache directory. [free text] """,
            [],
            fixedchoices=False,
            default=f"{str(pathlib.Path.home())}/.viroconstrictor/containers",
        )

    # subprocess.call("/bin/clear", shell=False)

    with open(file, "w") as conffile:
        conf_object.write(conffile)
        log.info("[green]Successfully written global configuration settings[/green]")


def AllOptionsGiven(config: configparser.ConfigParser) -> bool:
    """Function checks if all required config options are present in the already existing config file.
    Necessary to avoid missing config options when a user updates to a new version of ViroConstrictor.

    Parameters
    ----------
    config
        The configuration file.

    Returns
    -------
        A boolean value.

    """
    all_present: bool = True

    if config.has_section("COMPUTING") is True:
        if (
            config.has_option("COMPUTING", "compmode") is True
            and config["COMPUTING"]["compmode"] == "grid"
            and config.has_option("COMPUTING", "queuename") is False
            or config.has_option("COMPUTING", "compmode") is not True
        ):
            all_present = False
    else:
        all_present = False

    if config.has_section("GENERAL") is True:
        if (
            config.has_option("GENERAL", "auto_update") is True
            and config["GENERAL"]["auto_update"] == "no"
            and config.has_option("GENERAL", "ask_for_update") is False
            or config.has_option("GENERAL", "auto_update") is not True
        ):
            all_present = False
    else:
        all_present = False

    if config.has_section("REPRODUCTION") is True:
        if (
            config.has_option("REPRODUCTION", "repro_method") is True
            and config["REPRODUCTION"]["repro_method"] == "containers"
            and config.has_option("REPRODUCTION", "container_cache_path") is False
            or config.has_option("REPRODUCTION", "repro_method") is not True
        ):
            all_present = False
    else:
        all_present = False

    return all_present


def ReadConfig(file: pathlib.Path) -> configparser.ConfigParser:
    """ReadConfig() reads a config file, and if it doesn't exist, it creates it.
    -> If it does exist, but is empty, the configfile is recreated.
    -> If it exists and is populated, it reads it.
    -> If it exists and is populated, but not all necessary options are given, the configfile is recreated

    Parameters
    ----------
    file
        The file to read from.

    Returns
    -------
        A configparser object

    """
    if FileExists(file) is False:
        BuildConfig(file)
    if FileExists(file) is True and FileIsPopulated(file) is False:
        BuildConfig(file)

    config = configparser.ConfigParser()
    config.read(file)

    while AllOptionsGiven(config) is False:
        BuildConfig(file)
        config = configparser.ConfigParser()
        config.read(file)
    log.info("[green]Succesfully read global configuration file[/green]")
    return config
