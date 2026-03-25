"""Tests for :mod:`ViroConstrictor.functions`.

This module validates the custom argument formatter and tab-completion helpers
used by the command-line interface.
"""

from __future__ import annotations

from argparse import SUPPRESS
from types import SimpleNamespace

import pytest

import ViroConstrictor.functions as vc_functions


def test_flexible_arg_formatter_formats_help_without_crashing_on_small_terminal(monkeypatch: pytest.MonkeyPatch) -> None:
    """Verify FlexibleArgFormatter handles very narrow terminals without crashing.

    Tests that the formatter correctly wraps and displays help text even when
    the terminal width is extremely limited (10 columns).

    Parameters
    ----------
    monkeypatch : pytest.MonkeyPatch
        Pytest fixture for mocking module attributes.
    """
    with monkeypatch.context() as local_patch:
        local_patch.setattr(vc_functions.shutil, "get_terminal_size", lambda: SimpleNamespace(columns=10))
        parser = vc_functions.RichParser(formatter_class=vc_functions.FlexibleArgFormatter)
        parser.add_argument("--threads", default=4, help="Number of worker threads")

        rendered_help = parser.format_help()

    assert "--threads" in rendered_help
    assert "Number" in rendered_help
    assert "worker" in rendered_help
    assert "threads" in rendered_help


def test_flexible_arg_formatter_formats_help_without_crashing_on_wide_terminal(monkeypatch: pytest.MonkeyPatch) -> None:
    """Verify FlexibleArgFormatter handles very wide terminals without crashing.

    Tests that the formatter correctly handles extremely wide terminal widths
    (220 columns) without format errors.

    Parameters
    ----------
    monkeypatch : pytest.MonkeyPatch
        Pytest fixture for mocking module attributes.
    """
    with monkeypatch.context() as local_patch:
        local_patch.setattr(vc_functions.shutil, "get_terminal_size", lambda: SimpleNamespace(columns=220))
        parser = vc_functions.RichParser(formatter_class=vc_functions.FlexibleArgFormatter)
        parser.add_argument("--target", default="SARS", help="Virus target")

        rendered_help = parser.format_help()

    assert "--target" in rendered_help
    assert "Virus target" in rendered_help


def test_get_help_string_appends_default_value_when_not_present() -> None:
    """Verify _get_help_string appends default values to help text.

    Tests that the FlexibleArgFormatter adds "(default: <value>)" when the
    help text doesn't already document the default.
    """
    parser = vc_functions.RichParser(formatter_class=vc_functions.FlexibleArgFormatter)
    parser.add_argument("--threads", default=4, help="Number of worker threads")

    action = next(a for a in parser._actions if "--threads" in a.option_strings)
    formatter = vc_functions.FlexibleArgFormatter("prog")
    help_text = formatter._get_help_string(action)

    assert help_text is not None
    assert "Number of worker threads" in help_text
    assert "default: 4" in help_text


@pytest.mark.parametrize(
    "help_text,default,expected",
    [
        ("Default is already documented", 2, "Default is already documented"),
        ("value", SUPPRESS, "value"),
        ("value", None, "value"),
        (None, 3, None),
    ],
)
def test_get_help_string_does_not_append_when_not_needed(help_text: str | None, default: object, expected: str | None) -> None:
    """Verify _get_help_string skips appending defaults in appropriate cases.

    Tests that _get_help_string does not append defaults when:
    - The help text already documents the default
    - The default is SUPPRESS
    - The default is None
    - The help text itself is None

    Parameters
    ----------
    help_text : str | None
        The help text provided to an argument.
    default : object
        The default value for the argument.
    expected : str | None
        The expected help text returned by _get_help_string.
    """
    action = SimpleNamespace(help=help_text, default=default)
    formatter = vc_functions.FlexibleArgFormatter("prog")
    assert formatter._get_help_string(action) == expected


def test_format_action_invocation_optional_argument_uses_metavar_once() -> None:
    """Verify _format_action_invocation shows metavar only once for optional args.

    For optional arguments with short and long forms (e.g., -t, --threads),
    the metavar should appear only once in the formatted invocation string,
    not duplicated for each option form.
    """
    parser = vc_functions.RichParser(formatter_class=vc_functions.FlexibleArgFormatter)
    parser.add_argument("-t", "--threads", default=2)

    action = next(a for a in parser._actions if "--threads" in a.option_strings)
    formatter = vc_functions.FlexibleArgFormatter("prog")
    invocation = formatter._format_action_invocation(action)

    assert "-t" in invocation
    assert "--threads" in invocation
    assert invocation.count("THREADS") == 1


def test_format_action_invocation_falls_back_for_positional_argument() -> None:
    """Verify _format_action_invocation falls back gracefully for positional args.

    For positional arguments (non-optional), _format_action_invocation should
    return the argument name without formatting errors or extra punctuation.
    """
    parser = vc_functions.RichParser(formatter_class=vc_functions.FlexibleArgFormatter)
    parser.add_argument("sample")

    action = next(a for a in parser._actions if a.dest == "sample")
    formatter = vc_functions.FlexibleArgFormatter("prog")
    invocation = formatter._format_action_invocation(action)

    assert "sample" in invocation
    assert "," not in invocation


def test_format_action_invocation_falls_back_for_store_const_action() -> None:
    """Verify _format_action_invocation handles store_const actions correctly.

    For store_const actions (boolean flags), the formatter should return just
    the flag name (e.g., "--flag") without metavar or extra formatting.
    """
    parser = vc_functions.RichParser(formatter_class=vc_functions.FlexibleArgFormatter)
    parser.add_argument("--flag", action="store_const", const=True)

    action = next(a for a in parser._actions if "--flag" in a.option_strings)
    formatter = vc_functions.FlexibleArgFormatter("prog")
    invocation = formatter._format_action_invocation(action)

    assert invocation == "--flag"


def test_split_lines_and_fill_text_wrap_paragraphs_and_preserve_bullets() -> None:
    """Verify _split_lines and _fill_text preserve bullet points and wrap text.

    Tests that the FlexibleArgFormatter correctly:
    - Wraps long text to fit narrow widths
    - Preserves bullet point markers (*, -) at paragraph starts
    - Handles multi-line bullet items with proper continuation indentation
    """
    formatter = vc_functions.FlexibleArgFormatter("prog")
    text = """
        Intro line that should wrap into multiple pieces because the width is intentionally tiny.

        * first bullet item with more words than fit the target width
          continuation line for first bullet
    """

    split_lines = formatter._split_lines(text, width=35)
    filled_text = formatter._fill_text(text, width=35, indent="")

    assert len(split_lines) > 2
    assert split_lines == filled_text.split("\n")
    assert any(line.lstrip().startswith("*") for line in split_lines)


@pytest.mark.parametrize(
    "line,expected",
    [
        ("plain text", (0, 0)),
        ("    plain text", (4, 4)),
        ("  * bullet text", (2, 4)),
        ("  1) numbered", (2, 5)),
        ("\tstarts with tab", (0, 0)),
    ],
)
def test_indents_detects_regular_and_list_lines(line: str, expected: tuple[int, int]) -> None:
    """Verify _indents correctly identifies indent levels for various line types.

    Tests that _indents returns a tuple of (current_indent, continuation_indent)
    for plain text, bulleted lists (*, -), numbered lists (1), 2)), etc.), and
    tab-indented lines.

    Parameters
    ----------
    line : str
        A line of text to analyze for indentation.
    expected : tuple[int, int]
        Expected (current_indent, continuation_indent) tuple.
    """
    formatter = vc_functions.FlexibleArgFormatter("prog")
    assert formatter._indents(line) == expected


def test_split_paragraphs_merges_same_indent_lines_and_keeps_blank_lines() -> None:
    """Verify _split_paragraphs merges lines and preserves blank lines.

    Tests that _split_paragraphs:
    - Merges consecutive lines with the same indentation level
    - Preserves blank lines as paragraph separators
    - Maintains bullet list structure and indentation
    """
    formatter = vc_functions.FlexibleArgFormatter("prog")
    text = """
        alpha
        beta

          * one
            continuation
    """

    paragraphs = formatter._split_paragraphs(text)

    assert paragraphs[0] == "alpha beta"
    assert "" in paragraphs
    assert paragraphs[2].lstrip().startswith("*")


def test_para_reformat_returns_placeholder_for_blank_input() -> None:
    """Verify _para_reformat returns a list for blank/whitespace-only input.

    Tests that _para_reformat safely handles edge cases like empty strings or
    whitespace-only input without raising exceptions, always returning a list
    of strings.
    """
    formatter = vc_functions.FlexibleArgFormatter("prog")
    reformatted = formatter._para_reformat("\n\n", width=20)
    assert isinstance(reformatted, list)
    assert all(isinstance(x, str) for x in reformatted)


def test_rich_parser_print_message_uses_rich_print(monkeypatch: pytest.MonkeyPatch) -> None:
    """Verify RichParser._print_message delegates to rich.print.

    Tests that RichParser correctly integrates with the Rich library for
    formatted console output, using rich.print instead of standard argparse
    printing.

    Parameters
    ----------
    monkeypatch : pytest.MonkeyPatch
        Pytest fixture for mocking module attributes.
    """
    parser = vc_functions.RichParser()
    seen: list[str] = []
    monkeypatch.setattr(vc_functions.rich, "print", lambda message: seen.append(message))

    parser._print_message("hello", file=None)

    assert seen == ["hello"]


@pytest.mark.xfail(reason="Known tab completion behavior bug pending fix")
def test_path_completer_expands_home_and_preserves_suffix_before_globbing(monkeypatch: pytest.MonkeyPatch) -> None:
    """Verify pathCompleter expands ~ without dropping the typed suffix.

    Tests that tabCompleter.pathCompleter:
    - Expands a user path like ~/sam to /home/user/sam (prefix preserved)
    - Adds trailing slash to directory completions
    - Returns indexed matches from glob results

    Parameters
    ----------
    monkeypatch : pytest.MonkeyPatch
        Pytest fixture for mocking module attributes.
    """
    completer = vc_functions.tabCompleter()
    seen_patterns: list[str] = []

    monkeypatch.setattr(vc_functions.readline, "get_line_buffer", lambda: "~/sam")
    monkeypatch.setattr(vc_functions.os.path, "expanduser", lambda value: value.replace("~", "/tmp/home", 1))
    monkeypatch.setattr(vc_functions.os.path, "isdir", lambda path: path == "/tmp/home/sam")

    def fake_glob(pattern: str) -> list[str]:
        seen_patterns.append(pattern)
        return [f"{pattern}A", f"{pattern}B"]

    monkeypatch.setattr(vc_functions.glob, "glob", fake_glob)

    first = completer.pathCompleter("~/sam", 0)
    second = completer.pathCompleter("~/sam", 1)

    assert seen_patterns == ["/tmp/home/sam/*", "/tmp/home/sam/*"]
    assert first.startswith("/tmp/home/sam/")
    assert second.startswith("/tmp/home/sam/")
    assert first.endswith("A")
    assert second.endswith("B")


@pytest.mark.xfail(reason="Known tab completion behavior bug pending fix")
def test_path_completer_returns_none_for_out_of_range_state(monkeypatch: pytest.MonkeyPatch) -> None:
    """Verify pathCompleter returns None for out-of-range states.

    Tests that tabCompleter.pathCompleter:
    - Returns indexed glob match results
    - Returns None when index exceeds available matches (readline protocol)
    - Appends wildcard (*) to incomplete paths before globbing

    Parameters
    ----------
    monkeypatch : pytest.MonkeyPatch
        Pytest fixture for mocking module attributes.
    """
    completer = vc_functions.tabCompleter()
    monkeypatch.setattr(vc_functions.readline, "get_line_buffer", lambda: "abc")
    monkeypatch.setattr(vc_functions.os.path, "isdir", lambda _path: False)
    monkeypatch.setattr(vc_functions.glob, "glob", lambda pattern: [f"{pattern}1"])

    assert completer.pathCompleter("abc", 0) == "abc*1"
    assert completer.pathCompleter("abc", 1) is None


@pytest.mark.xfail(reason="Known tab completion behavior bug pending fix")
def test_create_list_completer_uses_line_prefix_when_line_present(monkeypatch: pytest.MonkeyPatch) -> None:
    """Verify listCompleter filters matches by current line prefix.

    Tests that tabCompleter.listCompleter:
    - Filters list items by the current line buffer prefix
    - Returns indexed matches in order matching the prefix
    - Returns None for out-of-range indices (readline protocol)

    Parameters
    ----------
    monkeypatch : pytest.MonkeyPatch
        Pytest fixture for mocking module attributes.
    """
    completer = vc_functions.tabCompleter()
    completer.createListCompleter(["alpha", "beta", "alpine"])
    monkeypatch.setattr(vc_functions.readline, "get_line_buffer", lambda: "al")

    assert completer.listCompleter("ignored", 0) == "alpha"
    assert completer.listCompleter("ignored", 1) == "alpine"
    assert completer.listCompleter("ignored", 2) is None


@pytest.mark.xfail(reason="Known tab completion behavior bug pending fix")
def test_create_list_completer_returns_options_with_trailing_space_for_empty_line(monkeypatch: pytest.MonkeyPatch) -> None:
    """Verify listCompleter appends trailing space to matches for empty buffer.

    Tests that tabCompleter.listCompleter:
    - Returns all list items when the line buffer is empty
    - Appends trailing space to each completion for readline integration
    - Returns None for out-of-range indices (readline protocol)

    Parameters
    ----------
    monkeypatch : pytest.MonkeyPatch
        Pytest fixture for mocking module attributes.
    """
    completer = vc_functions.tabCompleter()
    completer.createListCompleter(["alpha", "beta"])
    monkeypatch.setattr(vc_functions.readline, "get_line_buffer", lambda: "")

    assert completer.listCompleter("ignored", 0) == "alpha "
    assert completer.listCompleter("ignored", 1) == "beta "
    assert completer.listCompleter("ignored", 2) is None
