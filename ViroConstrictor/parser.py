import argparse
import multiprocessing
import os
import pathlib
import re
import sys
from typing import Any, Hashable, List

import numpy as np
import pandas as pd

from ViroConstrictor import __prog__, __version__
from ViroConstrictor.functions import FlexibleArgFormatter, RichParser
from ViroConstrictor.logging import log, setup_logger
from ViroConstrictor.samplesheet import GetSamples
from ViroConstrictor.userprofile import ReadConfig
from ViroConstrictor.validatefasta import CheckReferenceFile
from ViroConstrictor.workflow.presets import match_preset_name


class CLIparser:
    def __init__(self, input_args: list[str]) -> None:
        self.flags: argparse.Namespace = self._get_args(input_args)
        self.logfile = setup_logger(self.flags.output)
        log.info(f"ViroConstrictor version: [blue]{__version__}[/blue]")
        self.cli_errors = self._validate_cli_args()
        if self.cli_errors:
            for err in self.cli_errors:
                log.error(err)
            sys.exit(1)
        self.user_config = ReadConfig(
            pathlib.Path("~/.ViroConstrictor_defaultprofile.ini").expanduser()
        )
        self.flags.presets = self.flags.disable_presets is False
        self.samples_df = pd.DataFrame()
        self.samples_dict: dict[Hashable, Any] = {}
        if self.flags.samplesheet is not None:
            self._print_missing_asset_warning(self.flags, True)
            self.samples_dict = self._make_samples_dict(
                self._check_sample_sheet(self.flags.samplesheet),
                self.flags,
                GetSamples(self.flags.input, self.flags.platform),
            )
            self.samples_df = pd.DataFrame.from_dict(self.samples_dict, orient="index")
        else:
            self._print_missing_asset_warning(self.flags, False)
            self.samples_dict = self._make_samples_dict(
                None, self.flags, GetSamples(self.flags.input, self.flags.platform)
            )
        if self.samples_df.empty:
            self.samples_df = pd.DataFrame.from_dict(self.samples_dict, orient="index")
        self.samples_df = self.samples_df.reset_index(drop=False).rename(
            columns={"index": "SAMPLE"}
        )
        (
            self.input_path,
            self.workdir,
            self.exec_start_path,
            self.snakefile,
        ) = self._get_paths_for_workflow(self.flags)
        if not self.samples_dict:
            sys.exit(1)
        log.info("[green]Successfully parsed all command line arguments[/green]")
        self._check_sample_properties(
            self.samples_dict
        )  # raises errors if stuff is not right

    def _validate_cli_args(self) -> list[str] | None:
        arg_errors = []
        if dir_path(self.flags.input) is False:
            arg_errors.append(
                f"'[magenta]{self.flags.input}[/magenta]' is not a directory."
            )
        if self.flags.samplesheet is not None:
            allowed_extensions = [".xls", ".xlsx", ".csv", ".tsv"]
            if file_exists(self.flags.samplesheet) is False:
                arg_errors.append(
                    f"'[magenta]{self.flags.samplesheet}[/magenta]' is not an existing file."
                )
            if (
                check_file_extension(
                    allowed_extensions=allowed_extensions, fname=self.flags.samplesheet
                )
                is False
            ):
                arg_errors.append(
                    f"'[magenta]{self.flags.samplesheet}[/magenta]' does not have a valid file extension.\nAllowed file extenstions for the samplesheet: {' '.join([f'[blue]{x}[/blue]' for x in allowed_extensions])}"
                )
        if self.flags.reference is not None:
            allowed_extensions = [".fasta", ".fa"]
            if self.flags.reference == "NONE":
                arg_errors.append(
                    f"'[magenta]{self.flags.reference}[/magenta]' cannot be given for the reference file."
                )
            if file_exists(self.flags.reference) is False:
                arg_errors.append(
                    f"'[magenta]{self.flags.reference}[/magenta]' is not an existing file."
                )
            if (
                check_file_extension(
                    allowed_extensions=allowed_extensions, fname=self.flags.reference
                )
                is False
            ):
                arg_errors.append(
                    f"'[magenta]{self.flags.reference}[/magenta]' does not have a valid file extension.\nAllowed file extenstions for the reference: {' '.join([f'[blue]{x}[/blue]' for x in allowed_extensions])}"
                )
        if self.flags.primers is not None:
            allowed_extensions = [".fasta", ".fa", ".bed"]
            if file_exists(self.flags.primers) is False:
                arg_errors.append(
                    f"'[magenta]{self.flags.primers}[/magenta]' is not an existing file."
                )
            if (
                check_file_extension(
                    allowed_extensions=allowed_extensions, fname=self.flags.primers
                )
                is False
            ):
                arg_errors.append(
                    f"'[magenta]{self.flags.primers}[/magenta]' does not have a valid file extension.\nAllowed file extenstions for the primers: {' '.join([f'[blue]{x}[/blue]' for x in allowed_extensions])}"
                )
        if self.flags.features is not None:
            allowed_extensions = [".gff", ".gff3"]
            if file_exists(self.flags.features) is False:
                arg_errors.append(
                    f"'[magenta]{self.flags.features}[/magenta]' is not an existing file."
                )
            if (
                check_file_extension(
                    allowed_extensions=allowed_extensions, fname=self.flags.features
                )
                is False
            ):
                arg_errors.append(
                    f"'[magenta]{self.flags.features}[/magenta]' does not have a valid file extension.\nAllowed file extenstions for the features: {' '.join([f'[blue]{x}[/blue]' for x in allowed_extensions])}"
                )
        return arg_errors

    def _check_sample_properties(self, sampleinfo: dict[Hashable, Any]) -> None:
        """Check that the reference fasta file exists and is valid

        Parameters
        ----------
        sampleinfo : dict[str, dict[str, str]]
            A dictionary of dictionaries. The outer dictionary is keyed by sample name, and the inner
        dictionary is keyed by the parameter name.

        """
        reference_files: set = set()
        for item in sampleinfo:
            if sample := sampleinfo.get(item):
                if reffile := sample.get("REFERENCE"):
                    if not os.path.isfile(reffile):
                        log.error(
                            f"[bold red]The given reference fasta file for sample '{item}' does not exist. Please check the reference fasta and try again. Exiting...[/bold red]"
                        )
                        sys.exit(1)
                    reference_files.add(reffile)
        reference_files = set(reference_files)
        for f in reference_files:
            CheckReferenceFile(f)

    def _get_args(self, givenargs: list[str]) -> argparse.Namespace:
        """
        Parse the commandline args
        """
        parser: argparse.ArgumentParser = RichParser(
            prog=f"[bold]{__prog__}[/bold]",
            usage="%(prog)s \[required arguments] \[optional arguments]",
            description="%(prog)s: a pipeline for analysing Viral targeted (amplicon) sequencing data in order to generate a biologically valid consensus sequence.",
            formatter_class=FlexibleArgFormatter,
            add_help=False,
        )
        required_args = parser.add_argument_group("Required arguments")
        optional_args = parser.add_argument_group("Optional arguments")

        required_args.add_argument(
            "--input",
            "-i",
            type=str,
            metavar="DIR",
            help="The input directory with raw fastq(.gz) files",
            required=True,
        )

        required_args.add_argument(
            "--output",
            "-o",
            metavar="DIR",
            type=str,
            default=os.getcwd(),  # Default output dir is the current working dir
            help="Output directory",
            required=True,
        )

        optional_args.add_argument(
            "--samplesheet",
            "-samples",
            metavar="File",
            type=str,
            # type=lambda s: check_file_extension([".xls", ".xlsx", ".csv", ".tsv"], s),
            help="Sample sheet information file",
        )

        optional_args.add_argument(
            "--reference",
            "-ref",
            type=str,
            # type=lambda s: check_file_extension([".fasta", ".fa"], s),
            metavar="File",
            help="Input Reference sequence genome in FASTA format",
        )

        optional_args.add_argument(
            "--primers",
            "-pr",
            type=str,
            # type=lambda s: check_file_extension([".fasta", ".fa", ".bed"], s),
            metavar="File",
            help="Used primer sequences in FASTA or BED format. If no primers should be removed, supply the value NONE to this flag.",
        )

        required_args.add_argument(
            "--platform",
            default="nanopore",
            const="nanopore",
            nargs="?",
            choices=("nanopore", "illumina", "iontorrent"),
            help="Define the sequencing platform that was used to generate the dataset, either being 'nanopore', 'illumina' or 'iontorrent', see the docs for more info",
            required=True,
            metavar="'nanopore'/'illumina'/'iontorrent'",
        )

        required_args.add_argument(
            "--amplicon-type",
            "-at",
            default="end-to-end",
            const="end-to-end",
            nargs="?",
            choices=("end-to-end", "end-to-mid", "fragmented"),
            help="Define the amplicon-type, either being 'end-to-end', 'end-to-mid', or 'fragmented'. See the docs for more info",
            required=True,
            metavar="'end-to-end'/'end-to-mid'/'fragmented'",
        )

        optional_args.add_argument(
            "--target",
            "--preset",
            metavar="Str",
            help="Define the specific target for the pipeline, if the target matches a certain preset then pre-defined analysis settings will be used, see the docs for more info",
        )

        optional_args.add_argument(
            "--match-ref",
            "-mr",
            default=False,
            action="store_true",
            help="Match your data to the best reference available in the given reference fasta file.",
        )

        optional_args.add_argument(
            "--min-coverage",
            "-mc",
            default=30,
            type=int,
            metavar="N",
            help="Minimum coverage for the consensus sequence.",
        )

        optional_args.add_argument(
            "--features",
            "-gff",
            type=str,
            # type=lambda s: check_file_extension([".gff", ".gff3"], s),
            metavar="File",
            help="GFF file containing the Open Reading Frame (ORF) information of the reference. Supplying NONE will let ViroConstrictor use prodigal to determine coding regions",
        )

        optional_args.add_argument(
            "--primer-mismatch-rate",
            "-pmr",
            type=float,
            default=0.1,
            metavar="N",
            help="Maximum number of mismatches allowed in the primer sequences during primer coordinate search. Use 0 for exact primer matches\nDefault is 3.",
        )

        optional_args.add_argument(
            "--disable-presets",
            "-dp",
            action="store_true",
            default=False,
            help="Disable the use of presets, this will cause all analysis settings to be set to default values",
        )

        optional_args.add_argument(
            "--threads",
            "-t",
            default=min(multiprocessing.cpu_count(), 128),
            metavar="N",
            type=int,
            help=f"Number of local threads that are available to use.\nDefault is the number of available threads in your system ({min(multiprocessing.cpu_count(), 128)})",
        )

        optional_args.add_argument(
            "--version",
            "-v",
            version=__version__,
            action="version",
            help="Show the ViroConstrictor version and exit",
        )

        optional_args.add_argument(
            "--help",
            "-h",
            action="help",
            default=argparse.SUPPRESS,
            help="Show this help message and exit",
        )

        optional_args.add_argument(
            "--dryrun",
            action="store_true",
            help="Run the workflow without actually doing anything",
        )

        optional_args.add_argument(
            "--skip-updates",
            action="store_true",
            help="Skip the update check",
        )

        if not givenargs:
            log.error(
                f"{parser.prog} was called but no arguments were given, please try again\nUse '[cyan]{parser.prog} -h[/cyan]' to see the help document"
            )
            sys.exit(1)
        else:
            flags = parser.parse_args(givenargs)

        return flags

    def _check_sample_sheet(self, file: str) -> pd.DataFrame:
        """Checks the sample sheet for the required columns.

        Parameters
        ----------
        df
            the sample sheet

        Returns
        -------
            the result of the check_samplesheet_rows function. (pandas dataframe)

        """
        df = open_sample_sheet(file)
        if not df.empty:
            df.columns = df.columns.str.upper()
            req_cols = check_samplesheet_columns(df)
            if req_cols is False:
                sys.exit(1)
            df = samplesheet_enforce_absolute_paths(df)
            if df.get("PRESET") is None:
                df[["PRESET", "PRESET_SCORE"]] = df.apply(
                    lambda x: pd.Series(
                        match_preset_name(x["VIRUS"], use_presets=self.flags.presets)
                    ),
                    axis=1,
                )
            return check_samplesheet_rows(df)
        return pd.DataFrame()

    def _make_samples_dict(
        self,
        df: pd.DataFrame | None,
        args: argparse.Namespace,
        filedict: dict[str, str] | dict[str, dict[str, str]],
    ) -> dict[Hashable, Any]:
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
            log.error(
                f"'[magenta]{args.input}[/magenta]' does not contain any valid FastQ files. Exiting..."
            )
            sys.exit(1)
        log.info(
            f"[green]Valid FastQ files were found in the input directory.[/green] ('[magenta]{args.input}[/magenta]')"
        )
        indirFrame: pd.DataFrame = sampledir_to_df(filedict, args.platform)
        if df is not None:
            df.set_index("SAMPLE", inplace=True)
            df = pd.merge(df, indirFrame, left_index=True, right_index=True)
            if df.empty:
                log.error(
                    "[bold red]The files given in the samplesheet do not match the files given in the input-directory. Please check your samplesheet or input directory and try again.[/bold red]"
                )
                sys.exit(1)
            if len(indirFrame) > len(df):
                log.error(
                    "[bold red]Not all samples in the input directory are present in the given samplesheet. Please check your samplesheet or input directory and try again.[/bold red]"
                )
                sys.exit(1)
            if len(indirFrame) < len(df):
                log.error(
                    "[bold red]Either not all samples in the samplesheet are present in the given input directory, or there are duplicate samples in the samplesheet. Please check your samplesheet or input directory and try again.[/bold red]"
                )
                sys.exit(1)
            if df.get("PRIMER-MISMATCH-RATE") is None:
                df["PRIMER-MISMATCH-RATE"] = args.primer_mismatch_rate
            if df.get("MIN-COVERAGE") is None:
                df["MIN-COVERAGE"] = args.min_coverage
            if df.get("PRIMERS") is None:
                df["PRIMERS"] = args.primers
                if args.primers is None:
                    log.error(
                        "[bold red]No primer file specified in samplesheet or in command line options. Consider adding the -pr flag.[/bold red]"
                    )
                    sys.exit(1)
            if df.get("FEATURES") is None:
                df["FEATURES"] = args.features
                if args.features is None:
                    log.error(
                        "[bold red]No features file specified in samplesheet or in command line options. Consider adding the -gff flag.[/bold red]"
                    )
                    sys.exit(1)
            if df.get("PRESET") is None:
                df[["PRESET", "PRESET_SCORE"]] = df.apply(
                    lambda x: pd.Series(
                        match_preset_name(x["VIRUS"], use_presets=args.presets)
                    ),
                    axis=1,
                )
            df = pd.DataFrame.replace(df, np.nan, None)
            return df.to_dict(orient="index")
        return args_to_df(args, indirFrame).to_dict(orient="index")

    def _print_missing_asset_warning(
        self, args: argparse.Namespace, sheet_present: bool
    ) -> None:
        """If a sample sheet is present, print a warning that conflicting run-wide settings given through the commandline will be ignored.
        If no sample sheet is present, check if all required run-wide settings are given. If not, exit with a corresponding error message

        Parameters
        ----------
        args : argparse.Namespace
            argparse.Namespace
        sheet_present : bool
            boolean, whether a sample sheet was provided

        """
        if sheet_present:
            if args.primers is not None:
                log.warning(
                    "[yellow]Both a sample sheet and run-wide primer file was given, the primer file given through the commandline will be ignored[/yellow]"
                )
            if args.reference is not None:
                log.warning(
                    "[yellow]Both a sample sheet and run-wide reference fasta was given, the reference fasta given through the commandline will be ignored[/yellow]"
                )
            if args.features is not None:
                log.warn(
                    "[yellow]Both a sample sheet and run-wide GFF file was given, the GFF file given through the commandline will be ignored[/yellow]"
                )
        if not sheet_present and any(
            map(
                lambda f: f is None,
                {args.primers, args.reference, args.features, args.target},
            )
        ):
            log.error(
                f"[bold red]Run-wide analysis settings were not provided and no samplesheet was given either with per-sample run information.\nPlease either provide all required information ([underline]reference[/underline], [underline]primers[/underline], [underline]genomic features[/underline] and [underline]viral-target[/underline]) for a run-wide analysis or provide a samplesheet with per-sample run information[/bold red]"
            )
            sys.exit(1)

    def _get_paths_for_workflow(
        self, flags: argparse.Namespace
    ) -> tuple[str, str, str, str]:
        """Takes the input and output paths from the command line, and then creates the working directory if
        it doesn't exist. It then changes the current working directory to the working directory

        Parameters
        ----------
        flags : argparse.Namespace
            argparse.Namespace: The flags that were passed to the script.

        Returns
        -------
            A tuple of strings.

        """
        input_path: str = os.path.abspath(flags.input)
        working_directory: str = os.path.abspath(flags.output)
        exec_start_path: str = os.path.abspath(os.getcwd())
        snakefile: str = os.path.join(
            os.path.abspath(os.path.dirname(__file__)), "workflow", "workflow.smk"
        )

        if not os.path.exists(working_directory):
            os.makedirs(working_directory)
        if os.getcwd() != working_directory:
            os.chdir(working_directory)

        return input_path, working_directory, exec_start_path, snakefile


def samplesheet_enforce_absolute_paths(df: pd.DataFrame) -> pd.DataFrame:
    """
    Ensure the columns in the dataframe which may contain paths are always absolute paths.
    This is necessary as the pipeline may be run from different working directories.
    The columns are "PRIMERS", "FEATURES", "REFERENCE".
    If the value is "NONE", it is left as is.
    If the value is a relative path, it is converted to an absolute path.
    If the value is an absolute path, it is left as is.
    If there is a '~' in the path, it is expanded to the user's home directory.

    Parameters
    ----------
    df : pd.DataFrame
        The dataframe to enforce absolute paths on.

    Returns
    -------
    pd.DataFrame
        The modified dataframe with enforced absolute paths.

    Raises
    ------
    TypeError
        If the input argument is not a pandas DataFrame.

    """
    if not isinstance(df, pd.DataFrame):
        raise TypeError("Input argument must be a pandas DataFrame.")
    columns_to_enforce: List[str] = ["PRIMERS", "FEATURES", "REFERENCE"]
    for column in columns_to_enforce:
        if column in df.columns:
            df[column] = df[column].apply(
                lambda x: os.path.abspath(os.path.expanduser(x)) if x != "NONE" else x
            )
    return df


def file_exists(path: str) -> bool:
    """Check if a file exists.

    Parameters
    ----------
    path : str
        The path to the file.

    Returns
    -------
    bool
        A boolean value.

    """
    return True if path == "NONE" else os.path.isfile(path)


def is_excel_file(ext: str) -> bool:
    """If the extension is in the list of Excel extensions, return True, otherwise return False.

    Only checks the file extension, not file integrity.

    Parameters
    ----------
    ext : str
        The extension of the file.

    Returns
    -------
    bool
        A boolean value.

    """
    return ext in {".xls", ".xlsx"}


def is_csv_file(ext: str) -> bool:
    """Return True if the file extension is .csv, otherwise return False.

    Only checks the file extension, not file integrity

    Parameters
    ----------
    ext : str
        The extension of the file.

    Returns
    -------
    bool
        A boolean value.

    """
    return ext in {".csv"}


def is_tsv_file(ext: str) -> bool:
    """If the extension is in the list of extensions, return True, otherwise return False.

    Only checks the file extension, not file integrity.

    Parameters
    ----------
    ext : str
        the file extension

    Returns
    -------
    bool
        True if the extension is in the list of extensions, False otherwise
    """
    return ext in {".tsv"}


def open_sample_sheet(file: str) -> pd.DataFrame:
    """Given a file, return a pandas dataframe created with correct open function

    Parameters
    ----------
    file: Path
        The file to open.

    Returns
    -------
    pd.DataFrame
        A pandas dataframe.

    """
    # check if file is not empty
    if os.stat(file).st_size == 0:
        log.error("[red]Samplesheet file is empty.[/red]")
        sys.exit(1)
    file_extension = "".join(pathlib.Path(file).suffixes)
    try:
        if is_excel_file(file_extension):
            return pd.read_excel(file)
        if is_csv_file(file_extension):
            return pd.read_csv(file)
        if is_tsv_file(file_extension):
            return pd.read_csv(file, sep="\t")
    except Exception:
        log.exception(f"{file} is not a valid samplesheet file type.")
        return pd.DataFrame()
    return pd.DataFrame()


def required_cols(cols: list[str]) -> bool:
    """The function required_cols takes a list of column names and checks to see if the required columns
    are present

    Parameters
    ----------
    cols : list
        a list of column names

    Returns
    -------
    bool
        A boolean value.

    """
    cols = [c.upper() for c in cols]
    return all(i in cols for i in ["SAMPLE", "VIRUS", "REFERENCE"])


def check_samplesheet_columns(df: pd.DataFrame) -> bool:
    """Wrapper-function to check whether the samplesheet file has all the required columns or not

    Parameters
    ----------
    df
        the dataframe of the samplesheet

    Returns
    -------
    bool
        A boolean value.

    """
    if not required_cols(df.columns.tolist()):
        log.error("[bold red]Missing required columns in samplesheet file.[/bold red]")
        return False
    return True


def check_samplesheet_rows(df: pd.DataFrame) -> pd.DataFrame:
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
        "PRESET": {
            "dtype": str,
            "required": True,
            "disallowed_characters": None,
            "path": False,
        },
        "PRESET_SCORE": {
            "dtype": float,
            "required": True,
            "disallowed_characters": None,
            "path": False,
        },
    }
    for colName, colValue in df.items():
        if colName not in formats:
            log.error(
                f"[bold red]Unknown column '{colName}' in samplesheet.[/bold red]\n[yellow]Please check the column-headers in your samplesheet file and try again.\nAllowed column-headers are as follows: {' | '.join(list(formats))}[/yellow]"
            )
            sys.exit(1)
        if formats[colName]["required"] is True:
            if any(colValue.isnull()):
                log.error(
                    f"[bold red]Not all required information is given in column [underline]{colName}[/underline][/bold red]\n[yellow]Please check your samplesheet file and try again.[/yellow]"
                )
                sys.exit(1)
            for val in colValue:
                if not isinstance(val, formats[colName]["dtype"]):
                    log.error(
                        f"[bold red]{colName} column contains invalid data type.[/bold red]\n[yellow]Please check your samplesheet file and try again.[/yellow]"
                    )
                    sys.exit(1)
                if (
                    formats[colName]["disallowed_characters"] is not None
                    and formats[colName]["dtype"] == str
                ):
                    chars = re.compile(formats[colName]["disallowed_characters"])
                    if chars.search(val):
                        log.error(
                            f"[bold red]{colName} column contains one or more invalid characters.[/bold red]\n[yellow]Please check your samplesheet file and try again.[/yellow]"
                        )
                        sys.exit(1)
                if formats[colName]["path"] is True and not os.path.isfile(val):
                    log.error(
                        f"[bold red]{colName} column contains a path which doesn't point to a file: '{val}'.[/bold red]\n[yellow]Please check your samplesheet file and try again.[/yellow]"
                    )
                    sys.exit(1)
    return df


def check_file_extension(allowed_extensions: list[str], fname: str) -> bool:
    """If the input file name is "NONE", return it; otherwise, check that the file exists and has a valid
    extension, and return the absolute path to the file

    Parameters
    ----------
    fname
        The name of the file to be checked.
    allowed_extensions
        a list of file extensions that are allowed

    Returns
    -------
        The absolute path of the file.

    """
    if fname == "NONE":
        return True
    ext = "".join(pathlib.Path(fname).suffixes)
    if not any(ext.endswith(c) for c in allowed_extensions):
        return False
    return True


def dir_path(arginput: str) -> bool:
    """If the input is a directory, return it. Otherwise, print an error message and exit

    Parameters
    ----------
    arginput : str
        The input directory.

    Returns
    -------
    str
        the directory path.

    """
    return bool(os.path.isdir(arginput))


def CheckInputFiles(indir: str) -> bool:
    """Check if the input files are valid fastq files

    The function takes one argument, indir, which is the directory where the input files are located.
    The function returns a boolean value, True or False depending on whether the files in 'indir' are valid

    Parameters
    ----------
    indir : str
        The directory where the input files are located

    Returns
    -------
    bool
        A boolean value.

    """
    allowedextensions: list[str] = [".fastq", ".fq", ".fastq.gz", ".fq.gz"]
    foundfiles: list[str] = []

    for filenames in os.listdir(indir):
        extensions = "".join(pathlib.Path(filenames).suffixes)
        foundfiles.append(extensions)

    return any(
        file
        for file in foundfiles
        if any(file.endswith(ext) for ext in allowedextensions)
    )


def args_to_df(args: argparse.Namespace, df: pd.DataFrame) -> pd.DataFrame:
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
    df["PRIMERS"] = os.path.abspath(args.primers) if args.primers != "NONE" else "NONE"
    df["REFERENCE"] = os.path.abspath(args.reference)
    df["FEATURES"] = (
        os.path.abspath(args.features) if args.features != "NONE" else "NONE"
    )
    df["MIN-COVERAGE"] = args.min_coverage
    df["PRIMER-MISMATCH-RATE"] = args.primer_mismatch_rate
    df[["PRESET", "PRESET_SCORE"]] = match_preset_name(args.target, args.presets)
    df = pd.DataFrame.replace(df, np.nan, None)
    return df


def sampledir_to_df(
    sampledict: dict[str, str] | dict[str, dict[str, str]], platform: str
) -> pd.DataFrame:
    """Takes a dictionary of sample names and lists of input files, and returns a dataframe with the
    sample names as the index and the input files as the columns

    Parameters
    ----------
    sampledict : dict[str, list[str]]
        A dictionary of sample names to lists of input files.
    platform : str
        The sequencing platform used to generate the data.

    Returns
    -------
        A dataframe with the sample names as the index and the input files as the columns.

    """
    frame = pd.DataFrame.from_dict(sampledict, orient="index")
    if platform == "illumina":
        frame.index.rename("SAMPLE", inplace=True)
        return frame
    if platform in {"nanopore", "iontorrent"}:
        frame.index.rename("SAMPLE", inplace=True)
        frame.rename(columns={0: "INPUTFILE"}, inplace=True)
        return frame
    raise ValueError(f"Platform {platform} not supported")
