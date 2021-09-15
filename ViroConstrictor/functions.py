# pylint: disable=C0103.C0301

"""
Basic functions for various uses throughout ViroConstrictor
"""

import argparse
import glob
import os
import shutil

import readline


class MyHelpFormatter(argparse.RawTextHelpFormatter):
    """
    This is a custom formatter class for argparse. It allows for some custom formatting,
    in particular for the help texts with multiple options (like bridging mode and verbosity level).
    http://stackoverflow.com/questions/3853722
    """

    def __init__(self, prog):
        terminal_width = shutil.get_terminal_size().columns
        os.environ["COLUMNS"] = str(terminal_width)
        max_help_position = min(max(24, terminal_width // 2), 80)
        super().__init__(prog, max_help_position=max_help_position)

    def _get_help_string(self, action):
        help_text = action.help
        if (
            action.default != argparse.SUPPRESS
            and "default" not in help_text.lower()
            and action.default is not None
        ):
            help_text += " (default: " + str(action.default) + ")"
        return help_text


class color:
    """
    define basic colors to use in the terminal
    """

    PURPLE = "\033[95m"
    CYAN = "\033[96m"
    DARKCYAN = "\033[36m"
    BLUE = "\033[94m"
    GREEN = "\033[92m"
    YELLOW = "\033[93m"
    RED = "\033[91m"
    BOLD = "\033[1m"
    UNDERLINE = "\033[4m"
    END = "\033[0m"


# tabCompleter Class taken from https://gist.github.com/iamatypeofwalrus/5637895
## this was intended for the raw_input() function of python. But that one is deprecated now
## this also seems to work for the new input() functions however so muy bueno
#! use the gnureadline module instead of readline module
##! the 'normal' readline module causes a bug with memory pointers
##! --> https://stackoverflow.com/questions/43013060/python-3-6-1-crashed-after-i-installed-readline-module
class tabCompleter:
    """
    A tab completer that can either complete from
    the filesystem or from a list.

    Partially taken from:
    http://stackoverflow.com/questions/5637124/tab-completion-in-pythons-raw-input
    """

    def pathCompleter(self, text, state):
        """
        This is the tab completer for systems paths.
        Only tested on *nix systems
        """
        line = readline.get_line_buffer().split()

        # replace ~ with the user's home dir. See https://docs.python.org/2/library/os.path.html
        if "~" in text:
            text = os.path.expanduser("~")

        # autocomplete directories with having a trailing slash
        if os.path.isdir(text):
            text += "/"

        return [x for x in glob.glob(text + "*")][state]

    def createListCompleter(self, ll):
        """
        This is a closure that creates a method that autocompletes from
        the given list.

        Since the autocomplete function can't be given a list to complete from
        a closure is used to create the listCompleter function with a list to complete
        from.
        """

        def listCompleter(text, state):
            line = readline.get_line_buffer()

            if not line:
                return [c + " " for c in ll][state]

            return [c + " " for c in ll if c.startswith(line)][state]

        self.listCompleter = listCompleter
