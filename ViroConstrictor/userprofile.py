# pylint: disable=C0103

"""
Read or write the configuration file for ViroConstrictor.
This is so a user doesn't have to give this information manually on every run
"""
import configparser
import os
import readline
import subprocess
import sys

from .functions import color, tabCompleter


def FileExists(file):
    """Function returns a boolean value (True or False) depending on whether the file exists or not

    Parameters
    ----------
    file
        The file to check if it exists.

    Returns
    -------
        True or False

    """
    return bool(os.path.isfile(file))


def FileIsPopulated(file):
    """If the file exists and is not empty, return True. Otherwise, return False

    Parameters
    ----------
    file
        The file to check.

    Returns
    -------
        The size of the file in bytes.

    """
    return os.stat(file).st_size >= 1


def AskPrompts(intro, prompt, options, fixedchoices=False):
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
    if fixedchoices is True:
        completer = tabCompleter()
        completer.createListCompleter(options)

        readline.set_completer_delims("\t")
        readline.parse_and_bind("tab: complete")
        readline.set_completer(completer.listCompleter)

    subprocess.call("/bin/clear", shell=False)
    print(intro)
    while "the answer is invalid":
        if fixedchoices is True:
            reply = input(prompt).lower().strip()
            if reply in options:
                return reply
            if reply == "quit":
                print("Quitting...")
                sys.exit(-1)
            else:
                print(
                    "The given answer was invalid. Please choose one of the available options\n"
                )
        if fixedchoices is False:
            reply = input(prompt).strip()
            if reply == "quit":
                sys.exit(-1)
            else:
                return reply


def BuildConfig(file):
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

    conf_object["COMPUTING"] = {
        "compmode": AskPrompts(
            f"""
ViroConstrictor can run in two computing-modes.
{color.YELLOW + color.UNDERLINE}local{color.END} or {color.YELLOW + color.UNDERLINE}HPC/Grid{color.END}
Please specify the computing-mode that you wish to use for ViroConstrictor.
            """,
            f"""Do you wish to run ViroConstrictor in {color.YELLOW}local{color.END} or {color.YELLOW}grid{color.END} mode? [local/grid] """,
            ["local", "grid"],
            fixedchoices=True,
        )
    }

    if conf_object["COMPUTING"]["compmode"] == "grid":
        conf_object["COMPUTING"]["queuename"] = AskPrompts(
            f"""
Grid mode has been chosen. Please enter the name of computing-queue that you wish to use on your grid/HPC cluster.\nThis is necessary so ViroConstrictor will send all the various tasks to the correct (remote) computers.\n\n{color.BOLD + color.UNDERLINE + color.YELLOW}Please note that this is case-sensitive{color.END}\n""",
            "Please specify the name of the Queue on your grid/HPC cluster that you wish to use. [free text] ",
            [],
            fixedchoices=False,
        )

    conf_object["GENERAL"] = {
        "auto_update": AskPrompts(
            f"""
ViroConstrictor can check and update itself everytime you run it.
Please specify whether you wish to enable the auto-update feature.
            """,
            f"""Do you wish to enable the auto-update feature? [yes/no] """,
            ["yes", "no"],
            fixedchoices=True,
        )
    }

    if conf_object["GENERAL"]["auto_update"] == "no":
        conf_object["GENERAL"]["ask_for_update"] = AskPrompts(
            f"""
ViroConstrictor will not automatically update itself, but ViroConstrictor can still check for updates and ask you if you wish to update.
            """,
            f"""Do you want ViroConstrictor to {color.YELLOW}ask you{color.END} to update everytime a new update is available? [yes/no] """,
            ["yes", "no"],
            fixedchoices=True,
        )

    with open(file, "w") as conffile:
        conf_object.write(conffile)


def AllOptionsGiven(config):
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
    all_present = True

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

    return all_present


def ReadConfig(file):
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
    return config
