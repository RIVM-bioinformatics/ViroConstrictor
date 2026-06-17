"""Unit tests for BaseScript helper class.

These tests focus on intended script behavior: argument registration, command-line
entrypoint wiring, and failure modes.
"""

from __future__ import annotations

import sys
from argparse import ArgumentParser
from pathlib import Path

import pytest

from ViroConstrictor.workflow.helpers.base_script_class import BaseScript


class _RecordingScript(BaseScript):
    """Concrete test script that records constructor and run invocations."""

    constructed: dict[str, object] | None = None
    ran: bool = False

    @classmethod
    def reset(cls) -> None:
        """Reset shared state used by assertions."""
        cls.constructed = None
        cls.ran = False

    @classmethod
    def add_arguments(cls, parser: ArgumentParser) -> None:
        """Extend base arguments with an extra required option."""
        super().add_arguments(parser)
        parser.add_argument("--mode", required=True)

    def __init__(self, input: Path | str, output: Path | str, mode: str) -> None:
        """Initialize recording script and store constructor arguments for assertion.

        Parameters
        ----------
        input : Path | str
            Input file path to pass to parent BaseScript.
        output : Path | str
            Output file path to pass to parent BaseScript.
        mode : str
            Custom mode argument to store in shared class state for verification.
        """
        super().__init__(input=input, output=output)
        _RecordingScript.constructed = {
            "input": input,
            "output": output,
            "mode": mode,
        }

    def run(self) -> None:
        """Record that run was reached."""
        _RecordingScript.ran = True


class _CrashingScript(BaseScript):
    """Concrete test script whose run method fails intentionally."""

    def __init__(self, input: Path | str, output: Path | str) -> None:
        """Initialize crashing script for error propagation testing.

        Parameters
        ----------
        input : Path | str
            Input file path to pass to parent BaseScript.
        output : Path | str
            Output file path to pass to parent BaseScript.
        """
        super().__init__(input=input, output=output)

    def run(self) -> None:
        """Raise an error to verify propagation through main()."""
        raise RuntimeError("boom")


def test_init_stores_input_and_output_without_mutation(tmp_path: Path) -> None:
    """Verify constructor stores given input and output values as provided."""
    in_path = tmp_path / "in.fastq"
    out_path = tmp_path / "out.fastq"

    script = BaseScript(input=in_path, output=str(out_path))

    assert script.input == in_path
    assert script.output == str(out_path)


def test_run_on_base_class_requires_subclass_override() -> None:
    """Verify base run method enforces subclass implementation contract."""
    with pytest.raises(NotImplementedError, match="Subclasses should implement this method"):
        BaseScript("input.txt", "output.txt").run()


def test_add_arguments_registers_required_io_options() -> None:
    """Verify add_arguments creates required --input and --output options."""
    parser = ArgumentParser()
    BaseScript.add_arguments(parser)

    parsed = parser.parse_args(["--input", "a.txt", "--output", "b.txt"])

    assert parsed.input == "a.txt"
    assert parsed.output == "b.txt"


def test_add_arguments_rejects_missing_required_option() -> None:
    """Verify CLI parser exits when required I/O arguments are incomplete."""
    parser = ArgumentParser()
    BaseScript.add_arguments(parser)

    with pytest.raises(SystemExit):
        parser.parse_args(["--input", "a.txt"])


def test_main_constructs_subclass_with_parsed_args_and_runs(monkeypatch: pytest.MonkeyPatch) -> None:
    """Verify main() parses argv, instantiates subclass, and calls run()."""
    _RecordingScript.reset()
    monkeypatch.setattr(sys, "argv", ["prog", "--input", "in.txt", "--output", "out.txt", "--mode", "strict"])

    _RecordingScript.main()

    assert _RecordingScript.ran is True
    assert _RecordingScript.constructed == {
        "input": "in.txt",
        "output": "out.txt",
        "mode": "strict",
    }


def test_main_reports_parse_errors_for_missing_subclass_required_args(monkeypatch: pytest.MonkeyPatch) -> None:
    """Verify main() exits with parser error when subclass args are missing."""
    _RecordingScript.reset()
    monkeypatch.setattr(sys, "argv", ["prog", "--input", "in.txt", "--output", "out.txt"])

    with pytest.raises(SystemExit) as exc_info:
        _RecordingScript.main()

    assert exc_info.value.code == 2
    assert _RecordingScript.ran is False


def test_main_propagates_run_exception(monkeypatch: pytest.MonkeyPatch) -> None:
    """Verify exceptions in run() are not swallowed by main()."""
    monkeypatch.setattr(sys, "argv", ["prog", "--input", "in.txt", "--output", "out.txt"])

    with pytest.raises(RuntimeError, match="boom"):
        _CrashingScript.main()


@pytest.mark.xfail(
    strict=True,
    reason="Intended behavior from class documentation: constructor should reject missing input paths.",
)
def test_constructor_rejects_nonexistent_input_path(tmp_path: Path) -> None:
    """Intended behavior: nonexistent input should fail during script initialization."""
    missing_input = tmp_path / "missing.fastq"

    with pytest.raises((FileNotFoundError, ValueError)):
        BaseScript(input=missing_input, output=tmp_path / "out.fastq")


@pytest.mark.xfail(
    strict=True,
    reason="Intended behavior from class documentation: constructor should reject pre-existing output paths.",
)
def test_constructor_rejects_existing_output_path(tmp_path: Path) -> None:
    """Intended behavior: pre-existing output should be rejected to avoid silent overwrite."""
    existing_output = tmp_path / "result.fastq"
    existing_output.write_text("already-there", encoding="utf-8")

    with pytest.raises((FileExistsError, ValueError)):
        BaseScript(input=tmp_path / "in.fastq", output=existing_output)
