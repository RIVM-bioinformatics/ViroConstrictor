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
    if os.path.isfile(file):
        return True
    return False


def FileIsPopulated(file):
    if not os.stat(file).st_size >= 1:
        return False
    return True


def AskPrompts(intro, prompt, options, fixedchoices=False):
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
    all_present = True

    if config.has_section("COMPUTING") is True:
        if config.has_option("COMPUTING", "compmode") is True and config["COMPUTING"]["compmode"] == "grid" and config.has_option("COMPUTING", "queuename") is False or config.has_option("COMPUTING", "compmode") is not True:
            all_present = False
    else:
        all_present = False

    if config.has_section("GENERAL") is True:
        if config.has_option("GENERAL", "auto_update") is True and config["GENERAL"]["auto_update"] == "no" and config.has_option("GENERAL", "ask_for_update") is False or config.has_option("GENERAL", "auto_update") is not True:
            all_present = False
    else:
        all_present = False

    return all_present


def ReadConfig(file):
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
