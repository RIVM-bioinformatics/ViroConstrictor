# pylint: disable=C0103

"""
Read or write the configuration file for ViroConstrictor.
This is so a user doesn't have to give this information manually on every run
"""
import configparser
import os
import readline
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

    os.system("clear")
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
            f"""Grid mode has been chosen. Please enter the name of computing-queue that you wish to use on your grid/HPC cluster.\nThis is necessary so ViroConstrictor will send all the various tasks to the correct (remote) computers.\n\n{color.BOLD + color.UNDERLINE + color.YELLOW}Please note that this is case-sensitive{color.END}\n""",
            "Please specify the name of the Queue on your grid/HPC cluster that you wish to use. [free text] ",
            [],
            fixedchoices=False,
        )

    with open(file, "w") as conffile:
        conf_object.write(conffile)


def ReadConfig(file):
    if FileExists(file) is False:
        BuildConfig(file)
    if FileExists(file) is True and FileIsPopulated(file) is False:
        BuildConfig(file)

    config = configparser.ConfigParser()
    config.read(file)
    return config
