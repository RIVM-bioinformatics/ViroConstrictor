# pylint: disable=C0103.C0301

"""
Basic functions for various uses throughout ViroConstrictor
"""

import glob
import os
import re
import readline
import shutil
import textwrap
from argparse import SUPPRESS, Action, ArgumentParser, HelpFormatter
from typing import IO, Optional

import rich


class FlexibleArgFormatter(HelpFormatter):
    """
    A subclass of ArgumentParser.HelpFormatter that fixes spacing in the help text and respects bullet points.
    Especially useful for multi-line help texts combined with default values.

    This class has taken a lot of inspiration from the 'argparse-formatter' made by Dave Steele; https://github.com/davesteele/argparse_formatter

    Direct link to davesteele's original code that served as an inspiration for this class: https://github.com/davesteele/argparse_formatter/blob/a15d89a99e20b0cad4c389a2aa490c551cef4f9c/argparse_formatter/flexi_formatter.py

    ---
    This class helps to alleviate the following points of the ArgParse help formatting:
    * The help text will be aligned with the argument name/flags, instead of printing the help description on a newline
    * Adjusting the width of the help text in relationship to the width of the terminal to make sure there is enough space between the argument name/flags and the help text (thus not overloading the end-user with an unreadable wall of text)
    * Adding a default value to the help text (on a newline, and indented) if one is provided in the ArgParse constructor
    * Respecting bullet points in the help description
    * Respecting newlines in the help description (you may have to add a space after the newline to make sure it is properly catched by the formatter)
    * Respecting indentation in the help description (up to a certain degree)
    * Changes the behaviour of the metavar to be only printed once per long AND shorthand argument, instead of printing the metavar multiple times for every possible flag.
    """

    def __init__(self, prog: str) -> None:
        term_width = shutil.get_terminal_size().columns
        max_help_position = min(max(24, term_width // 2), 80)
        super().__init__(prog, max_help_position=max_help_position)

    def _get_help_string(self, action: Action) -> Optional[str]:
        """ """
        help_text: Optional[str] = action.help
        if (
            help_text is not None
            and action.default != SUPPRESS
            and "default" not in help_text.lower()
            and action.default is not None
        ):
            help_text += f"\n  ([underline]default: {str(action.default)}[/underline])"
        return help_text

    def _format_action_invocation(self, action: Action) -> str:
        """ """
        if not action.option_strings or action.nargs == 0:
            return super()._format_action_invocation(action)
        default = self._get_default_metavar_for_optional(action)
        args_string = self._format_args(action, default)
        return ", ".join(action.option_strings) + " " + args_string

    def _split_lines(self, text: str, width: int) -> list[str]:
        return self._para_reformat(text, width)

    def _fill_text(self, text: str, width: int, indent: str) -> str:
        lines = self._para_reformat(text, width)
        return "\n".join(lines)

    def _indents(self, line: str) -> tuple[int, int]:
        """Return line indent level and "sub_indent" for bullet list text."""

        matches = re.match(r"( *)", line)
        indent = len(matches[1]) if matches else 0
        if list_match := re.match(r"( *)(([*\-+>]+|\w+\)|\w+\.) +)", line):
            sub_indent = indent + len(list_match[2])
        else:
            sub_indent = indent
        return (indent, sub_indent)

    def _split_paragraphs(self, text: str) -> list[str]:
        """Split text in to paragraphs of like-indented lines."""

        text = textwrap.dedent(text).strip()
        text = re.sub("\n\n\n+", "\n\n", text)

        last_sub_indent: Optional[int] = None
        paragraphs: list[str] = []
        for line in text.splitlines():
            (indent, sub_indent) = self._indents(line)
            is_text = re.search(r"[^\s]", line) is not None

            if is_text and indent == sub_indent == last_sub_indent:
                paragraphs[-1] += f" {line}"
            else:
                paragraphs.append(line)

            last_sub_indent = sub_indent if is_text else None
        return paragraphs

    def _para_reformat(self, text: str, width: int) -> list[str]:
        """Reformat text, by paragraph."""

        paragraphs: list[str] = []
        for paragraph in self._split_paragraphs(text):
            (indent, sub_indent) = self._indents(paragraph)

            paragraph = self._whitespace_matcher.sub(" ", paragraph).strip()
            new_paragraphs: list[str] = textwrap.wrap(
                text=paragraph,
                width=width,
                initial_indent=" " * indent,
                subsequent_indent=" " * sub_indent,
            )

            # Blank lines get eaten by textwrap, put it back with [' ']
            paragraphs.extend(new_paragraphs or [" "])

        return paragraphs


class RichParser(ArgumentParser):
    """
    A subclass of `argparse.ArgumentParser` that overrides the `_print_message` method to use
    `rich.print` instead of `print`
    """

    def _print_message(self, message: str, file: Optional[IO[str]] = None) -> None:
        return rich.print(message)


# tabCompleter Class taken from https://gist.github.com/iamatypeofwalrus/5637895
## this was intended for the raw_input() function of python. But that one is deprecated now
## However, this also seems to work for the new input() functions
class tabCompleter:
    """
    A tab completer that can either complete from
    the filesystem or from a list.

    Partially taken from:
    http://stackoverflow.com/questions/5637124/tab-completion-in-pythons-raw-input
    """

    def pathCompleter(self, text: str, state: int) -> str:
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

        # we explicitly to a list comprehension here instead of a call to the constructor as the this would otherwise break the autocompletion functionality of paths.
        return [x for x in glob.glob(f"{text}*")][state]

    def createListCompleter(self, ll: list[str]) -> None:
        """
        This is a closure that creates a method that autocompletes from
        the given list.

        Since the autocomplete function can't be given a list to complete from
        a closure is used to create the listCompleter function with a list to complete
        from.
        """

        def listCompleter(text: str, state: int) -> Optional[str]:
            line: str = readline.get_line_buffer()

            return (
                [c for c in ll if c.startswith(line)][state]
                if line
                else [f"{c} " for c in ll][state]
            )

        self.listCompleter = listCompleter
