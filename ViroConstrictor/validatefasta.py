# pylint: disable=C0103

"""
Basic functions to see if a fasta is valid
"""

import re

from Bio import SeqIO

from ViroConstrictor.logging import log


def ContainsSpecials(seq: str) -> bool:
    """It takes a string as input and returns True if the string contains any characters other than the 20
    amino acids, '-', or '*'

    Parameters
    ----------
    seq
        the sequence to be checked

    Returns
    -------
    bool
        True or False

    """
    chars = re.compile("[^actgumrwsykvhdbnACTGUMRWSYKVHDBN-]")
    return chars.search(seq) is not None


def ContainsAmbiguities(seq: str) -> bool:
    """If the sequence contains any of the characters in the string "umrwsykvhdbnUMRWSYKVHDBN" (all possible nucleotide ambiguities), then return
    True, otherwise return False

    Parameters
    ----------
    seq : str
        the sequence to be checked

    Returns
    -------
        A boolean value.

    """
    chars = re.compile("[umrwsykvhdbnUMRWSYKVHDBN]")
    return chars.search(seq) is not None


def IsValidRef(inputfile: str) -> bool:
    """If the input file is a valid FASTA file, and none of the sequences in the file contain ambiguous
    characters, then the file is a valid reference file

    Parameters
    ----------
    inputfile : str
        The path to the reference file.

    Returns
    -------
    bool
        A boolean value.

    """
    if IsValidFasta(inputfile):
        return not any(
            ContainsAmbiguities(str(record.seq))
            for record in SeqIO.parse(inputfile, "fasta")
        )
    return False


def IsValidFasta(inputfile: str) -> bool:
    """Function takes a fasta file (path) as input and returns True if all sequences in the file are valid, and False if
    any sequence in given fasta is invalid

    Parameters
    ----------
    inputfile : str
        The input file to check.

    Returns
    -------
    bool
        A boolean value.

    """
    if inputfile == "NONE":
        return True
    results = [
        ContainsSpecials(str(record.seq)) for record in SeqIO.parse(inputfile, "fasta")
    ]

    return not any(results)


def CheckReferenceFile(referencefile: str, warnings_as_errors: bool = False) -> None:
    """Checks that the reference file is in the correct format, and that it doesn't contain any
    ambiguities

    Parameters
    ----------
    referencefile : str
        The path to the reference file.
    warnings_as_errors : bool, optional
        If True, warnings will be treated as errors.

    """
    errors: list[Exception] = []
    warnings: list[str] = []
    for record in SeqIO.parse(referencefile, "fasta"):
        try:
            check_ref_header(record.id)
        except Exception as e:
            errors.append(e)

        # Check whether there are stretches of ambiguities
        matches = re.split("[ACTGactg]", str(record.seq))
        if longer_than_four := [m for m in matches if len(m) > 4]:
            errors.append(
                ValueError(
                    f"In file {referencefile}, record {record.id} has stretches of ambiguities:\n"
                    f"\t{longer_than_four}"
                )
            )

        # Check whether there are any ambiguous nucleotides
        if ambiguities := sum(map(len, matches)):
            unique_ambiguities = "".join(set("".join(matches)))
            w = f"""[cyan]{ambiguities}[/cyan] Ambiguous nucleotides found in file [magenta]{referencefile}[/magenta] in record [blue]{record.id}[/blue]:\t[bold yellow]{unique_ambiguities}[/bold yellow]\nPlease check whether this is intended."""
            if warnings_as_errors:
                errors.append(Exception(w))
            else:
                warnings.append(w)
    if warnings:
        for w in warnings:
            log.warning(f"{w}")
    if errors:
        for e in errors:
            log.error(f"{e}")
        exit(1)


def check_ref_header(header: str) -> str | None:
    """Checks the header of a reference fasta file to make sure there are no blacklisted or forbidden characters.
    Returns the header if it is valid, otherwise exits the program with an error message.

    Parameters
    ----------
    header : str
        The header of the reference fasta file.

    Returns
    -------
    str
        the fasta header.

    """
    if not header:
        log.error("Reference fasta does not have a header. Please add it.")
        exit(1)
    blacklisted_characters = {"\\", "/", ":", "*", "?", '"', "<", ">", "|", "\0"}
    if found_in_blacklist := {c for c in header if c in blacklisted_characters}:
        log.error(
            f"Reference fasta header '[bold red]{header}[/bold red]' contains the following invalid characters [bold red]{found_in_blacklist}[/bold red]\nPlease change the fasta header for this reference."
        )
        exit(1)
    reserved_words_on_windows = {
        "CON",
        "PRN",
        "AUX",
        "NUL",
        "COM1",
        "COM2",
        "COM3",
        "COM4",
        "COM5",
        "COM6",
        "COM7",
        "COM8",
        "COM9",
        "LPT1",
        "LPT2",
        "LPT3",
        "LPT4",
        "LPT5",
        "LPT6",
        "LPT7",
        "LPT8",
        "LPT9",
    }
    if header in reserved_words_on_windows:
        log.error(
            f"Reference fasta header '[bold red]{header}[/bold red]' is a reserved word on the windows operating system.\nPlease change it as it may cause problems."
        )
        exit(1)
    if all(c == "." for c in header):
        log.error(
            f"Reference fasta header '[bold red]{header}[/bold red]' is not a valid fasta header.\nPlease change the fasta header for this reference and try again."
        )
        exit(1)
    return header
