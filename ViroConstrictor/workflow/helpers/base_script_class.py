import logging
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

    def __init__(self, input: Path | str, output: Path | str, log_level: str = "INFO") -> None:
        self.input = input
        self.output = output
        self.logger = self._setup_logging(log_level)

    def run(self) -> None:
        raise NotImplementedError("Subclasses should implement this method.")

    @classmethod
    def add_arguments(cls, parser: ArgumentParser) -> None:
        parser.add_argument(
            "--input",
            metavar="File",
            help="Input file path.",
            type=str,
            required=True,
        )
        parser.add_argument(
            "--output",
            metavar="File",
            help="Output file path.",
            type=str,
            required=True,
        )

        parser.add_argument(
            "--log-level",
            metavar="Level",
            help="Logging level (e.g., DEBUG, INFO, WARNING, ERROR, CRITICAL).",
            type=str,
            default="INFO",
        )

    @classmethod
    def main(cls) -> None:
        parser = ArgumentParser(description=f"Run the {cls.__name__} script.")
        cls.add_arguments(parser)
        args = parser.parse_args()

        # Pass parsed arguments to the script
        script = cls(**vars(args))
        script.run()

    def _setup_logging(self, log_level: str) -> logging.Logger:
        """
        Besides snakemake rules, there should be logging for each script.
        The snakemake logs are written to logs/{rule_name}.{sample_name}.log
        Scripts will be logged to logs/scripts/{script_name}.log
        """
        numeric_level = getattr(logging, log_level.upper(), None)
        if not isinstance(numeric_level, int):
            # logging is not set yet, so we use print here
            print(f"Warning: Invalid log level: {log_level}. Using INFO instead.")
            numeric_level = logging.INFO

        logger = logging.getLogger(self.__class__.__name__)
        logger.setLevel(numeric_level)

        log_dir = Path("logs/scripts")
        log_dir.mkdir(parents=True, exist_ok=True)
        log_file = log_dir / f"{self.__class__.__name__}.log"

        file_handler = logging.FileHandler(log_file, mode="a")  # a means append mode
        console_handler = logging.StreamHandler()

        formatter = logging.Formatter("[%(asctime)s] %(levelname)s in %(module)s: %(message)s")
        file_handler.setFormatter(formatter)
        console_handler.setFormatter(formatter)

        logger.addHandler(file_handler)
        logger.addHandler(console_handler)

        return logger

    def log(self, message: str, level: int = logging.INFO) -> None:
        """
        Log a message at the specified logging level.
        Default level is INFO.
        """
        self.logger.log(level, message)
