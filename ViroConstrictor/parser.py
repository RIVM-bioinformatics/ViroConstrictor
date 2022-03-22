import argparse
import multiprocessing
import os
import pathlib
import re
import sys

import pandas as pd

from ViroConstrictor import __prog__, __version__
from ViroConstrictor.functions import MyHelpFormatter, color
from ViroConstrictor.samplesheet import GetSamples, nanopore_sheet


def is_excel_file(ext):
    return ext in [".xls", ".xlsx"]


def is_csv_file(ext):
    return ext in [".csv"]


def is_tsv_file(ext):
    return ext in [".tsv"]


def is_json_file(ext):
    return ext in [".json"]


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
    if is_json_file(file_extension):
        return pd.read_json(file, orient="records")
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
    if any(
        i not in cols for i in ["SAMPLE", "VIRUS", "PRIMERS", "REFERENCE", "FEATURES"]
    ):
        return False
    return True


def check_samplesheet_columns(df):
    if not required_cols(df.columns):
        print(
            f"{color.RED + color.BOLD}Missing required columns in samplesheet file.{color.END}",
            file=sys.stderr,
        )
        return False
    return True


def check_samplesheet_rows(df):
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
            "path": True,
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
            "path": True,
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
    df = open_sample_sheet(file)
    df.columns = df.columns.str.upper()
    req_cols = check_samplesheet_columns(df)
    if req_cols is False:
        sys.exit(1)
    return check_samplesheet_rows(df)


def is_valid_samplesheet_file(f):
    if os.path.isfile(f):
        if "".join(pathlib.Path(f).suffixes) in [
            ".xls",
            ".xlsx",
            ".csv",
            ".tsv",
            ".json",
        ]:
            return os.path.abspath(f)
        raise argparse.ArgumentTypeError(f"{f} is not a valid samplesheet file type.")
    print("Sample sheet file not found")
    sys.exit(1)


def check_input(choices, fname):
    if fname == "NONE":
        return fname
    if os.path.isfile(fname):
        ext = "".join(pathlib.Path(fname).suffixes)
        if ext not in choices:
            raise argparse.ArgumentTypeError(
                f"Input file doesn't end with one of {choices}"
            )
        return fname
    print(f'"{fname}" is not a file. Exiting...')
    sys.exit(-1)


def dir_path(arginput):
    if os.path.isdir(arginput):
        return arginput
    print(f'"{arginput}" is not a directory. Exiting...')
    sys.exit(1)


def currentpath():
    return os.getcwd()

def CheckInputFiles(indir):
    """
    Check if the input files are valid fastq files
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
        type=lambda s: check_input((".fasta", ".fa"), s),
        metavar="File",
        help="Used primer sequences in FASTA format",
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
        required=True
    )

    parser.add_argument(
        "--target",
        "--preset",
        metavar="Str",
        help="Define the specific target for the pipeline, if the target matches a certain preset then pre-defined analysis settings will be used, see the docs for more info",
    )
    
    parser.add_argument(
        '--match-ref',
        '-mr',
        default=False,
        action='store_true',
        help='Match your data to the best reference available in the given reference fasta file.'
    )
    
    parser.add_argument(
        '--min-coverage',
        '-mc',
        default=30,
        type=int,
        metavar="N",
        help='Minimum coverage for the consensus sequence.'
    )
    
    parser.add_argument(
        "--features",
        "-gff",
        type=lambda s: check_input((".gff"), s),
        metavar="File",
        help="GFF file containing the Open Reading Frame (ORF) information of the reference",
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
        "--skip-updates", action="store_true", help="Skip the update check",
    )

    if len(givenargs) < 1:
        print(
            f"{parser.prog} was called but no arguments were given, please try again\n\tUse '{parser.prog} -h' to see the help document"
        )
        sys.exit(1)
    else:
        flags = parser.parse_args(givenargs)

    return flags

def args_to_dict(args, df):
    df['VIRUS'] = args.target
    df['MATCH-REF'] = args.match_ref
    df['PRIMERS'] = args.primers
    df['REFERENCE'] = args.reference
    df['FEATURES'] = args.features
    df['MIN-COVERAGE'] = args.min_coverage
    df['PRIMER-MISMATCH-RATE'] = args.primer_mismatch_rate
    return df

def sampledir_to_df(sampledict, platform):
    """
    Convert the samplesheet to a pandas dataframe
    """
    frame = pd.DataFrame.from_dict(sampledict, orient="index")
    if platform == "illumina":
        frame.index.rename("SAMPLE", inplace=True)
        frame.rename(columns={"R1": "INPUTFILE_R1", "R2": "INPUTFILE_R2"}, inplace=True)
        return frame
    if platform in ['nanopore', 'iontorrent']:
        frame.index.rename("SAMPLE", inplace=True)
        frame.rename(columns={0: "INPUTFILE"}, inplace=True)
        return frame
    raise ValueError(f"Platform {platform} not supported")


def make_sampleinfo_dict(df, args, filedict):
    if not CheckInputFiles(args.input):
        print(f"\n{color.RED + color.BOLD}'{args.input}' does not contain any valid FastQ files. Exiting...{color.END}\n")
        sys.exit(1)
    print(f"\n{color.GREEN}Valid input files were found in the input directory.{color.END} ('{args.input}')\n")
    indirFrame = sampledir_to_df(filedict, args.platform)
    if df is not None:
        df.set_index("SAMPLE", inplace=True)
        df = pd.merge(df, indirFrame, left_index=True, right_index=True)
        if df.empty:
            print(f"\n{color.RED + color.BOLD}The files given in the samplesheet do not match the files given in the input-directory. Please check your samplesheet or input directory and try again.{color.END}\n")
            sys.exit(1)
        if len(indirFrame) > len(df):
            print(f"\n{color.RED + color.BOLD}Not all samples in the input directory are present in the given samplesheet. Please check your samplesheet or input directory and try again.{color.END}\n")
            sys.exit(1)
        if len(indirFrame) < len(df):
            print(f"\n{color.RED + color.BOLD}Not all sample in the samplesheet are present in the given input directory. Please check your samplesheet or input directory and try again.{color.END}\n")
            sys.exit(1)
        return df.to_dict(orient="index")
    return args_to_dict(args, indirFrame).to_dict(orient="index")
    


def ValidArgs(sysargs):
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
            print(f"{color.YELLOW}Both a sample sheet and run-wide primer file was given, the given run-wide primer file will be ignored{color.END}")
            args.primers = None
        if args.reference is not None:
            print(f"{color.YELLOW}Both a sample sheet and run-wide reference fasta was given, the given run-wide reference fasta will be ignored{color.END}")
            args.reference = None
        if args.features is not None:
            print(f"{color.YELLOW}Both a sample sheet and run-wide GFF file was given, the given run-wide GFF file will be ignored{color.END}")
            args.features = None

        df = check_sample_sheet(args.samplesheet)
        sampleinfo = make_sampleinfo_dict(df, args, GetSamples(args.input, args.platform))
    else:
        if args.primers is None or args.reference is None or args.features is None or args.target is None:
            print(f"{color.RED + color.BOLD}Run-wide analysis settings were not provided and no samplesheet was given either with per-sample run information.\nPlease either provide all required information (reference, primers, genomic features and viral-target) for a run-wide analysis or provide a samplesheet with per-sample run information{color.END}")
            sys.exit(1)
        sampleinfo = make_sampleinfo_dict(None, args, GetSamples(args.input, args.platform))

    #print(sampleinfo)
    if not sampleinfo:
        print("wut")
        sys.exit(1)
    return args, sampleinfo
