from argparse import ArgumentParser
from pathlib import Path


class BaseScript:
    """
    A base class for creating standardized scripts with argument parsing and execution logic.

    This class provides a common structure for scripts, including argument parsing,
    input/output path validation, and a `run` method that must be implemented by subclasses.
    Subclasses can extend the argument parser by overriding the `add_arguments` method.

    Attributes
    ----------
    input_path : Path
        The path to the input file.
    output_path : Path
        The path to the output file.

    Methods
    -------
    _validate_paths()
        Validates the existence of the input path and ensures the output path does not already exist.
    run()
        Abstract method to be implemented by subclasses with the script's main logic.
    add_arguments(parser)
        Class method that adds arguments to the argument parser. Subclasses can override this to add custom arguments.
    main()
        Entry point for the script. Parses arguments, initializes the script, and calls the `run` method.

    Notes
    -----
    To use this class:
    1. Inherit from `BaseScript` in your script.
    2. Implement the `run` method with the script's main logic.
    3. Override the `add_arguments` class method if additional arguments are needed.
    4. Call `MyScript.main()` in the script's main block to execute it.

    Example
    -------
    Here's an example of how to create a subclass with additional arguments:

    >>> from pathlib import Path
    >>> from base_script_class import BaseScript
    >>>
    >>> class MyScript(BaseScript):
    >>>     def __init__(self, input_path: Path, output_path: Path, extra_arg: str) -> None:
    >>>         super().__init__(input_path, output_path)
    >>>         self.extra_arg = extra_arg
    >>>
    >>>     def run(self) -> None:
    >>>         print(f"Processing {self.input_path} with extra_arg={self.extra_arg}")
    >>>
    >>>     @classmethod
    >>>     def add_arguments(cls, parser):
    >>>         super().add_arguments(parser)
    >>>         parser.add_argument("--extra_arg", type=str, required=True, help="An extra argument.")
    >>>
    >>> if __name__ == "__main__":
    >>>     MyScript.main()
    """

    def __init__(self, input_path: Path, output_path: Path) -> None:
        self.input_path = input_path
        self.output_path = output_path
        self._validate_paths()

    def _validate_paths(self) -> None:
        if not self.input_path.exists():
            raise FileNotFoundError(f"Input path {self.input_path} does not exist.")
        if self.output_path.exists():
            raise FileExistsError(f"Output path {self.output_path} already exists.")

    def run(self) -> None:
        raise NotImplementedError("Subclasses should implement this method.")

    @classmethod
    def add_arguments(cls, parser: ArgumentParser) -> None:
        parser.add_argument(
            "--input",
            metavar="File",
            help="Input file path.",
            type=Path,
            required=True,
        )
        parser.add_argument(
            "--output",
            metavar="File",
            help="Output file path.",
            type=Path,
            required=True,
        )

    @classmethod
    def main(cls) -> None:
        parser = ArgumentParser(description=f"Run the {cls.__name__} script.")
        cls.add_arguments(parser)
        args = parser.parse_args()

        # Pass parsed arguments to the script
        script = cls(**vars(args))
        script.run()
