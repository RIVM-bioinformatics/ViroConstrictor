import argparse
import multiprocessing
import os
import pathlib
import re
import sys

import pandas as pd

from ViroConstrictor import __prog__, __version__
from ViroConstrictor.functions import MyHelpFormatter, color
from ViroConstrictor.samplesheet import GetSamples


def is_excel_file(ext):
    """If the extension is in the list of Excel extensions, return True, otherwise return False

    Only checks the file extension, not file integrity.

    Parameters
    ----------
    ext
        The extension of the file.

    Returns
    -------
        A boolean value.

    """
    return ext in [".xls", ".xlsx"]


def is_csv_file(ext):
    """Return True if the file extension is .csv, otherwise return False.

    Only checks the file extension, not file integrity

    Parameters
    ----------
    ext
        The extension of the file.

    Returns
    -------
        A boolean value.

    """
    return ext in [".csv"]


def is_tsv_file(ext):
    """If the extension is in the list of extensions, return True, otherwise return False.

    Only checks the file extension, not file integrity.

    Parameters
    ----------
    ext
        the file extension

    Returns
    -------
        A list of all the files in the directory that end with the extension .tsv

    """
    return ext in [".tsv"]


def open_sample_sheet(file):
    """Given a file, return a pandas dataframe created with correct open function

    Parameters
    ----------
    file
        The file to open.

    Returns
    -------
        A pandas dataframe.

    """
    file_extension = "".join(pathlib.Path(file).suffixes)
    if is_excel_file(file_extension):
        return pd.read_excel(file)
    if is_csv_file(file_extension):
        return pd.read_csv(file)
    if is_tsv_file(file_extension):
        return pd.read_csv(file, sep="\t")
    raise TypeError(f"{file} is not a valid samplesheet file type.")


def required_cols(cols):
    """The function required_cols takes a list of column names and checks to see if the required columns
    are present

    Parameters
    ----------
    cols
        a list of column names

    Returns
    -------
        A boolean value.

    """
    cols = [c.upper() for c in cols]
    if any(i not in cols for i in ["SAMPLE", "VIRUS", "REFERENCE"]):
        return False
    return True


def check_samplesheet_columns(df):
    """Wrapper-function to check whether the samplesheet file has all the required columns or not

    Parameters
    ----------
    df
        the dataframe of the samplesheet

    Returns
    -------
        A boolean value.

    """
    if not required_cols(df.columns):
        print(
            f"{color.RED + color.BOLD}Missing required columns in samplesheet file.{color.END}",
            file=sys.stderr,
        )
        return False
    return True


def check_samplesheet_rows(df):
    """Checks whether the row-based contents of the samplesheet dataframe are valid.

    Parameters
    ----------
    df
        The dataframe containing the samplesheet

    Returns
    -------
        A dataframe with the columns: SAMPLE, VIRUS, PRIMERS, REFERENCE, FEATURES, MATCH-REF, MIN-COVERAGE,
    PRIMER-MISMATCH-RATE

    """
    formats = {
        "SAMPLE": {
            "dtype": str,
            "required": True,
            "disallowed_characters": "[ @!#$%^&*()+=|\]\[{}':;,/?<>~`]",
            "path": False,
        },
        "VIRUS": {
            "dtype": str,
            "required": True,
            "disallowed_characters": "[ @!#$%^&*()+=|\]\[{}':;,/?<>~`]",
            "path": False,
        },
        "PRIMERS": {
            "dtype": str,
            "required": True,
            "disallowed_characters": None,
            "path": False,
        },
        "REFERENCE": {
            "dtype": str,
            "required": True,
            "disallowed_characters": None,
            "path": True,
        },
        "FEATURES": {
            "dtype": str,
            "required": True,
            "disallowed_characters": None,
            "path": False,
        },
        "MATCH-REF": {
            "dtype": bool,
            "required": True,
            "disallowed_characters": None,
            "path": False,
        },
        "MIN-COVERAGE": {
            "dtype": int,
            "required": False,
            "disallowed_characters": None,
            "path": False,
        },
        "PRIMER-MISMATCH-RATE": {
            "dtype": int,
            "required": False,
            "disallowed_characters": None,
            "path": False,
        },
    }
    for (colName, colValue) in df.items():
        if colName not in formats:
            print(
                f"{color.RED + color.BOLD}\tUnknown column '{colName}' in samplesheet.{color.END}\n\tPlease check the column-headers in your samplesheet file and try again.\n\tAllowed column-headers are as follows: {' | '.join(list(formats))}",
                file=sys.stderr,
            )
            sys.exit(1)
        if formats[colName]["required"] is True:
            if any(colValue.isnull()):
                print(
                    f"{color.RED + color.BOLD}\tNot all required information is given in column {color.UNDERLINE +  colName}{color.END}\n\tPlease check your samplesheet file and try again.",
                    file=sys.stderr,
                )
                sys.exit(1)
            for val in colValue:
                if not isinstance(val, formats[colName]["dtype"]):
                    print(
                        f"{color.RED + color.BOLD}\t{colName} column contains invalid data type.{color.END}\n\tPlease check your samplesheet file and try again.",
                        file=sys.stderr,
                    )
                    sys.exit(1)
                if (
                    formats[colName]["disallowed_characters"] is not None
                    and formats[colName]["dtype"] == str
                ):
                    chars = re.compile(formats[colName]["disallowed_characters"])
                    if chars.search(val):
                        print(
                            f"{color.RED + color.BOLD}\t{colName} column contains one or more invalid characters.{color.END}\n\tPlease check your samplesheet file and try again.",
                            file=sys.stderr,
                        )
                        sys.exit(1)
                if formats[colName]["path"] is True and not os.path.isfile(val):
                    print(
                        f"{color.RED + color.BOLD}\t{colName} column contains a path which doesn't point to a file: '{val}'.{color.END}\n\tPlease check your samplesheet file and try again.",
                        file=sys.stderr,
                    )
                    sys.exit(1)
    return df


def check_sample_sheet(file):
    """Wrapper function that takes a samplesheet and triggers the checks for required columns and required (row-based) values.

    Parameters
    ----------
    file
        the path to the sample sheet

    Returns
    -------
        the result of the check_samplesheet_rows function. (pandas dataframe)

    """
    df = open_sample_sheet(file)
    df.columns = df.columns.str.upper()
    req_cols = check_samplesheet_columns(df)
    if req_cols is False:
        sys.exit(1)
    return check_samplesheet_rows(df)


def is_valid_samplesheet_file(f):
    """If the file exists and has a valid extension, return the absolute path to the file

    Parameters
    ----------
    f
        the file path to the samplesheet

    Returns
    -------
        The absolute path of the file.

    """
    if os.path.isfile(f):
        if "".join(pathlib.Path(f).suffixes) in {
            ".xls",
            ".xlsx",
            ".csv",
            ".tsv",
            ".json",
        }:
            return os.path.abspath(f)
        raise argparse.ArgumentTypeError(f"{f} is not a valid samplesheet file type.")
    print("Sample sheet file not found")
    sys.exit(1)


def check_input(choices, fname):
    """If the input file name is "NONE", return it; otherwise, check that the file exists and has a valid
    extension, and return the absolute path to the file

    Parameters
    ----------
    choices
        a list of file extensions that are allowed
    fname
        The name of the file to be checked.

    Returns
    -------
        The absolute path of the file.

    """
    if fname == "NONE":
        return fname
    if os.path.isfile(fname):
        ext = "".join(pathlib.Path(fname).suffixes)
        if ext not in choices:
            raise argparse.ArgumentTypeError(
                f"Input file doesn't end with one of {choices}"
            )
        return os.path.abspath(fname)
    print(f'"{fname}" is not a file. Exiting...')
    sys.exit(-1)


def dir_path(arginput):
    """If the input is a directory, return it. Otherwise, print an error message and exit

    Parameters
    ----------
    arginput
        The input directory.

    Returns
    -------
        the directory path.

    """
    if os.path.isdir(arginput):
        return arginput
    print(f'"{arginput}" is not a directory. Exiting...')
    sys.exit(1)


def currentpath():
    """Returns the current working directory

    Returns
    -------
        The current working directory.

    """
    return os.getcwd()


def CheckInputFiles(indir):
    """Check if the input files are valid fastq files

    The function takes one argument, indir, which is the directory where the input files are located.
    The function returns a boolean value, True or False

    Parameters
    ----------
    indir
        The directory where the input files are located

    Returns
    -------
        A boolean value.

    """
    allowedextensions = [".fastq", ".fq", ".fastq.gz", ".fq.gz"]
    foundfiles = []

    for filenames in os.listdir(indir):
        extensions = "".join(pathlib.Path(filenames).suffixes)
        foundfiles.append(extensions)

    return bool(any(i in allowedextensions for i in foundfiles))


def get_args(givenargs, parser):
    """
    Parse the commandline args
    """

    parser.add_argument(
        "--input",
        "-i",
        type=dir_path,
        metavar="DIR",
        help="The input directory with raw fastq(.gz) files",
        required=True,
    )

    parser.add_argument(
        "--output",
        "-o",
        metavar="DIR",
        type=str,
        default=currentpath(),
        help="Output directory",
        required=True,
    )

    parser.add_argument(
        "--samplesheet",
        "-samples",
        metavar="File",
        type=is_valid_samplesheet_file,
        help="Sample sheet information file",
    )

    parser.add_argument(
        "--reference",
        "-ref",
        type=lambda s: check_input((".fasta", ".fa"), s),
        metavar="File",
        help="Input Reference sequence genome in FASTA format",
    )

    parser.add_argument(
        "--primers",
        "-pr",
        type=lambda s: check_input((".fasta", ".fa", ".bed"), s),
        metavar="File",
        help="Used primer sequences in FASTA or BED format. If no primers should be removed, supply the value NONE to this flag.",
    )

    parser.add_argument(
        "--platform",
        default="nanopore",
        const="nanopore",
        nargs="?",
        choices=("nanopore", "illumina", "iontorrent"),
        help="Define the sequencing platform that was used to generate the dataset, either being 'nanopore', 'illumina' or 'iontorrent', see the docs for more info",
        required=True,
    )

    parser.add_argument(
        "--amplicon-type",
        "-at",
        default="end-to-end",
        const="end-to-end",
        nargs="?",
        choices=("end-to-end", "end-to-mid"),
        help="Define the amplicon-type, either being 'end-to-end' or 'end-to-mid', see the docs for more info",
        required=True,
    )

    parser.add_argument(
        "--target",
        "--preset",
        metavar="Str",
        help="Define the specific target for the pipeline, if the target matches a certain preset then pre-defined analysis settings will be used, see the docs for more info",
    )

    parser.add_argument(
        "--match-ref",
        "-mr",
        default=False,
        action="store_true",
        help="Match your data to the best reference available in the given reference fasta file.",
    )

    parser.add_argument(
        "--min-coverage",
        "-mc",
        default=30,
        type=int,
        metavar="N",
        help="Minimum coverage for the consensus sequence.",
    )

    parser.add_argument(
        "--features",
        "-gff",
        type=lambda s: check_input((".gff"), s),
        metavar="File",
        help="GFF file containing the Open Reading Frame (ORF) information of the reference. Supplying NONE will let ViroConstrictor use prodigal to determine coding regions",
    )

    parser.add_argument(
        "--primer-mismatch-rate",
        "-pmr",
        type=float,
        default=0.1,
        metavar="N",
        help="Maximum number of mismatches allowed in the primer sequences during primer coordinate search. Use 0 for exact primer matches\nDefault is 3.",
    )

    parser.add_argument(
        "--threads",
        "-t",
        default=min(multiprocessing.cpu_count(), 128),
        metavar="N",
        type=int,
        help=f"Number of local threads that are available to use.\nDefault is the number of available threads in your system ({min(multiprocessing.cpu_count(), 128)})",
    )

    parser.add_argument(
        "--version",
        "-v",
        version=__version__,
        action="version",
        help="Show the ViroConstrictor version and exit",
    )

    parser.add_argument(
        "--help",
        "-h",
        action="help",
        default=argparse.SUPPRESS,
        help="Show this help message and exit",
    )

    parser.add_argument(
        "--dryrun",
        action="store_true",
        help="Run the workflow without actually doing anything",
    )

    parser.add_argument(
        "--skip-updates",
        action="store_true",
        help="Skip the update check",
    )

    if len(givenargs) < 1:
        print(
            f"{parser.prog} was called but no arguments were given, please try again\n\tUse '{parser.prog} -h' to see the help document"
        )
        sys.exit(1)
    else:
        flags = parser.parse_args(givenargs)

    return flags


def args_to_df(args, df):
    """It takes the arguments from the command line and places them into a dataframe

    Parameters
    ----------
    args
        the arguments from the command line
    df
        the dataframe that will be used to store the results

    Returns
    -------
        A dataframe with the arguments as columns.

    """
    df["VIRUS"] = args.target
    df["MATCH-REF"] = args.match_ref
    df["PRIMERS"] = args.primers
    df["REFERENCE"] = args.reference
    df["FEATURES"] = args.features
    df["MIN-COVERAGE"] = args.min_coverage
    df["PRIMER-MISMATCH-RATE"] = args.primer_mismatch_rate
    return df


def sampledir_to_df(sampledict, platform):
    """Function converts a dictionary of samples to a pandas dataframe

    Parameters
    ----------
    sampledict
        a dictionary of sample names and their input files
    platform
        The sequencing platform used to generate the data.

    Returns
    -------
        A dataframe with the sample name as the index and the input file as the column.

    """
    frame = pd.DataFrame.from_dict(sampledict, orient="index")
    if platform == "illumina":
        frame.index.rename("SAMPLE", inplace=True)
        return frame
    if platform in ["nanopore", "iontorrent"]:
        frame.index.rename("SAMPLE", inplace=True)
        frame.rename(columns={0: "INPUTFILE"}, inplace=True)
        return frame
    raise ValueError(f"Platform {platform} not supported")


def make_sampleinfo_dict(df, args, filedict):
    """It takes a samplesheet (dataframe) and a dictionary fastq files, and returns a dictionary of sample
    information

    Parameters
    ----------
    df
        the samplesheet dataframe
    args
        The arguments given to the script.
    filedict
        a dictionary of the files in the input directory

    Returns
    -------
        A dictionary of the samplesheet and input directory.

    """
    if not CheckInputFiles(args.input):
        print(
            f"\n{color.RED + color.BOLD}'{args.input}' does not contain any valid FastQ files. Exiting...{color.END}\n"
        )
        sys.exit(1)
    print(
        f"\n{color.GREEN}Valid input files were found in the input directory.{color.END} ('{args.input}')\n"
    )
    indirFrame = sampledir_to_df(filedict, args.platform)
    if df is not None:
        df.set_index("SAMPLE", inplace=True)
        df = pd.merge(df, indirFrame, left_index=True, right_index=True)
        if df.empty:
            print(
                f"\n{color.RED + color.BOLD}The files given in the samplesheet do not match the files given in the input-directory. Please check your samplesheet or input directory and try again.{color.END}\n"
            )
            sys.exit(1)
        if len(indirFrame) > len(df):
            print(
                f"\n{color.RED + color.BOLD}Not all samples in the input directory are present in the given samplesheet. Please check your samplesheet or input directory and try again.{color.END}\n"
            )
            sys.exit(1)
        if len(indirFrame) < len(df):
            print(
                f"\n{color.RED + color.BOLD}Either not all samples in the samplesheet are present in the given input directory, or there are duplicate samples in the samplesheet. Please check your samplesheet or input directory and try again.{color.END}\n"
            )
            sys.exit(1)
        if df.get("PRIMER-MISMATCH-RATE") is None:
            df["PRIMER-MISMATCH-RATE"] = args.primer_mismatch_rate
        if df.get("MIN-COVERAGE") is None:
            df["MIN-COVERAGE"] = args.min_coverage
        if df.get("PRIMERS") is None:
            df["PRIMERS"] = args.primers
            if args.primers is None:
                print(
                    f"\n{color.RED + color.BOLD}No primer file specified in samplesheet or in command line options. Consider adding the -pr flag.{color.END}\n"
                )
                sys.exit(1)
        if df.get("FEATURES") is None:
            df["FEATURES"] = args.features
            if args.features is None:
                print(
                    f"\n{color.RED + color.BOLD}No features file specified in samplesheet or in command line options. Consider adding the -gff flag.{color.END}\n"
                )
                sys.exit(1)

        return df.to_dict(orient="index")
    return args_to_df(args, indirFrame).to_dict(orient="index")


def ValidArgs(sysargs):
    """Wrapper function which takes the command line arguments and returns a dictionary with all the information needed to run the pipeline

    Parameters
    ----------
    sysargs
        the command line arguments

    Returns
    -------
        args, sampleinfo

    """
    parser = argparse.ArgumentParser(
        prog=__prog__,
        usage=f"{__prog__} [required options] [optional arguments]",
        description="ViroConstrictor: a pipeline for analysing Viral targeted (amplicon) sequencing data in order to generate a biologically valid consensus sequence.",
        formatter_class=MyHelpFormatter,
        add_help=False,
    )
    args = get_args(sysargs, parser)

    if args.samplesheet is not None:
        if args.primers is not None:
            print(
                f"{color.YELLOW}Both a sample sheet and run-wide primer file was given, the given run-wide primer file will be ignored{color.END}"
            )
        if args.reference is not None:
            print(
                f"{color.YELLOW}Both a sample sheet and run-wide reference fasta was given, the given run-wide reference fasta will be ignored{color.END}"
            )
        if args.features is not None:
            print(
                f"{color.YELLOW}Both a sample sheet and run-wide GFF file was given, the given run-wide GFF file will be ignored{color.END}"
            )

        df = check_sample_sheet(args.samplesheet)
        sampleinfo = make_sampleinfo_dict(
            df, args, GetSamples(args.input, args.platform)
        )
    else:
        if any(
            map(
                lambda f: f is None,
                {args.primers, args.reference, args.features, args.target},
            )
        ):
            print(
                f"{color.RED + color.BOLD}Run-wide analysis settings were not provided and no samplesheet was given either with per-sample run information.\nPlease either provide all required information (reference, primers, genomic features and viral-target) for a run-wide analysis or provide a samplesheet with per-sample run information{color.END}"
            )
            sys.exit(1)
        sampleinfo = make_sampleinfo_dict(
            None, args, GetSamples(args.input, args.platform)
        )

    if not sampleinfo:
        sys.exit(1)
    return args, sampleinfo
