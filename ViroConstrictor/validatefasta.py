# pylint: disable=C0103

"""
Basic functions to see if a fasta is valid
"""

import re

from Bio import SeqIO

from ViroConstrictor.functions import color


def ContainsSpecials(seq):

    chars = re.compile("[^actgumrwsykvhdbnACTGUMRWSYKVHDBN-]")

    if chars.search(seq) is None:
        return False
    return True


def ContainsAmbiguities(seq):
    chars = re.compile("[umrwsykvhdbnUMRWSYKVHDBN]")
    if chars.search(seq) is None:
        return False
    return True


def IsValidRef(inputfile):
    if IsValidFasta(inputfile):
        return not any(
            ContainsAmbiguities(str(record.seq))
            for record in SeqIO.parse(inputfile, "fasta")
        )
    return False


def IsValidFasta(inputfile):
    if inputfile == "NONE":
        return True
    results = [
        ContainsSpecials(str(record.seq)) for record in SeqIO.parse(inputfile, "fasta")
    ]

    if any(results):
        return False
    return True


def CheckReferenceFile(referencefile, warnings_as_errors=False):
    errors = []
    warnings = []
    for record in SeqIO.parse(referencefile, "fasta"):
        try:
            check_ref_header(record.id)
        except Exception as e:
            errors.append(e)

        matches = re.split("[ACTGactg]", str(record.seq))
        if longer_than_four := [m for m in matches if len(m) > 4]:
            errors.append(
                ValueError(
                    f"In file {referencefile}, record {record.id} has stretches of ambiguities:\n"
                    f"\t{longer_than_four}"
                )
            )

        if ambiguities := sum(map(len, matches)):
            unique_ambiguities = "".join(set("".join(matches)))
            w = Warning(
                f"{ambiguities} Ambiguous nucleotides found in file {referencefile} in record {record.id}:\n"
                f"\t{unique_ambiguities}\n"
                "Check whether this is intended."
            )
            if warnings_as_errors:
                errors.append(w)
            else:
                warnings.append(w)
    if warnings:
        print(warnings)
    if warnings:
        print(color.YELLOW)
        for w in warnings:
            print(w.args[0])
        print(color.END)
    if errors:
        print(color.RED)
        print("Error:")
        for e in errors:
            print(e.args[0])
        print(color.END)
        exit(1)


def check_ref_header(s):
    if not s:
        raise ValueError("Reference fasta does not have a header. Please add it.")
    blacklisted_characters = {"\\", "/", ":", "*", "?", '"', "<", ">", "|", "\0"}
    if found_in_blacklist := [c for c in s if c in blacklisted_characters]:
        raise ValueError(
            f"Reference fasta header\n\t{s}\ncontains invalid characters\n\t{found_in_blacklist}\nPlease change the fasta header for this reference."
        )
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
    if s in reserved_words_on_windows:
        raise ValueError(
            f"Reference fasta header\n\t{s}\nis a reserved word on the windows operating system. Please change it."
        )
    if all(c == "." for c in s):
        raise ValueError(
            f"Reference fasta header\n\t{s}\nis not valid.\nPlease change the fasta header for this reference."
        )
    return s
