"""
This module provides functions for parsing and validating command-line arguments for the amplicon_covs script.

The script requires the following input files:
- A BED file with primers as given by AmpliGone.
- A TSV file with coverages as given by TrueConsense.
- A sample ID.
- An output TSV file for average coverage per amplicon.

The module includes the following functions:
- _common_file_checks: Perform common checks on a file to ensure it exists, is not empty, and is readable.
- check_primers_file: Validate the primers file to ensure it meets the required criteria.
- check_coverages_file: Validate the coverages file to ensure it meets the required criteria.
- check_output_file: Validate the output file to ensure it meets the required criteria.
- parse_args: Parse and validate command-line arguments for the script.

Each function raises an ArgumentTypeError if the input does not meet the required criteria.

Example usage:
--------------
To use this module, import it and call the parse_args function with the appropriate arguments.

    from amplicon_arg_parser import parse_args

    args = parse_args()
    print(args.primers)
    print(args.coverages)
    print(args.key)
    print(args.output)
"""

import os
from argparse import ArgumentParser, ArgumentTypeError, Namespace


def _common_file_checks(filename: str) -> None:
    """
    Perform common checks on a file to ensure it exists, is not empty, and is readable.

    This function raises an ArgumentTypeError if any of the following conditions are not met:
    - The file exists.
    - The file is not empty.
    - The file is readable.

    Parameters
    ----------
    filename : str
        The path to the file to be checked.

    Raises
    ------
    ArgumentTypeError
        If the file does not exist, is empty, or is not readable.
    """
    if not os.path.isfile(filename):
        raise ArgumentTypeError(f"File '{filename}' does not exist.")

    if os.path.getsize(filename) == 0:
        raise ArgumentTypeError(f"File '{filename}' is empty.")

    if not os.access(filename, os.R_OK):
        raise ArgumentTypeError(f"File '{filename}' is not readable.")


def check_primers_file(filename: str) -> str:
    """
    Validate the primers file to ensure it meets the required criteria.
    Does not open the file or check its contents.

    This function performs the following checks:
    - The file exists, is not empty, and is readable (using _common_file_checks).
    - The file has a .bed extension.

    Parameters
    ----------
    filename : str
        The path to the primers file to be checked.

    Returns
    -------
    str
        The validated filename.

    Raises
    ------
    ArgumentTypeError
        If the file does not exist, is empty, is not readable, or does not have a .bed extension.
    """
    _common_file_checks(filename)
    if not filename.lower().endswith(".bed"):
        raise ArgumentTypeError("Primers file must be a BED file.")

    return filename


def check_coverages_file(filename: str) -> str:
    """
    Validate the coverages file to ensure it meets the required criteria.
    Does not open the file or check its contents.

    This function performs the following checks:
    - The file exists, is not empty, and is readable (using _common_file_checks).
    - The file has a .tsv extension.

    Parameters
    ----------
    filename : str
        The path to the coverages file to be checked.

    Returns
    -------
    str
        The validated filename.

    Raises
    ------
    ArgumentTypeError
        If the file does not exist, is empty, is not readable, or does not have a .tsv extension.
    """
    _common_file_checks(filename)

    if not filename.lower().endswith(".tsv"):
        raise ArgumentTypeError("Coverages file must be a TSV file.")

    return filename


def check_output_file(filename: str) -> str:
    """
    Validate the output file to ensure it meets the required criteria.
    Does not open the file or check its contents.

    This function performs the following checks:
    - The file has a .csv extension.
    - The file does not already exist.
    - The directory containing the file is writable.

    Parameters
    ----------
    filename : str
        The path to the output file to be checked.

    Returns
    -------
    str
        The validated filename.

    Raises
    ------
    ArgumentTypeError
        If the file does not have a .csv extension, already exists, or the directory is not writable.
    """
    if not filename.lower().endswith(".csv"):
        raise ArgumentTypeError("Output file must be a CSV file.")

    if os.path.isfile(filename):
        raise ArgumentTypeError(
            f"Output file '{filename}' already exists. Please choose another name."
        )

    if not os.access(os.path.dirname(filename), os.W_OK):
        raise ArgumentTypeError(
            f"Directory '{os.path.dirname(filename)}' is not writable."
        )

    return filename


def parse_args(args: list[str] | None = None) -> Namespace:
    """
    Parse command-line arguments for the amplicov_covs script.

    This function sets up the argument parser and defines the required arguments for the script.
    It validates the input files using the specified check functions.

    Parameters
    ----------
    args : list[str], optional
        A list of command-line arguments to parse. If None, the arguments are taken from sys.argv.
        This parameter is used for testing purposes.

    Returns
    -------
    Namespace
        An argparse.Namespace object containing the parsed arguments.

    Arguments
    ---------
    --primers : File
        Input BED file with primers as given by AmpliGone. This file is validated by check_primers_file.

    --coverages : File
        Input file with coverages as given by TrueConsense. This file is validated by check_coverages_file.

    --key : String
        Sample ID.

    --output : File
        Output file with average coverage per amplicon. This file is validated by check_output_file.

    Raises
    ------
    ArgumentTypeError
        If any of the input files do not meet the required criteria.
    """
    parser = ArgumentParser()

    parser.add_argument(
        "--primers",
        metavar="File",
        type=check_primers_file,
        help="input BED file with primers as given by AmpliGone",
        required=True,
    )

    parser.add_argument(
        "--coverages",
        metavar="File",
        type=check_coverages_file,
        help="Input file with coverages as given by TrueConsense",
        required=True,
    )

    parser.add_argument(
        "--key",
        metavar="String",
        type=str,
        help="Sample ID",
        required=True,
    )

    parser.add_argument(
        "--output",
        metavar="File",
        type=check_output_file,
        help="Output file with average coverage per amplicon",
        required=True,
    )

    return parser.parse_args(args)
