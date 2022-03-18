import argparse
import re

from itertools import product

import pandas as pd
from Bio import SeqIO
from Bio import Seq


def FindAmbigousOptions(seq):
    ambigs = Seq.IUPACData.ambiguous_dna_values
    return list(map("".join, product(*map(ambigs.get, seq))))


def get_coords(seq, ref_seq, err_rate=0.1):
    max_errors = int(len(seq) * err_rate)

    options = FindAmbigousOptions(seq.upper())

    matches = set()
    for e in range(max_errors + 1):
        for option in options:
            for match in re.finditer(
                f"(?:{str(option)}){{s<={e}}}", ref_seq, re.IGNORECASE, concurrent=True
            ):
                matches.add(match.span())
    if matches:
        return matches
    print(f"Primer not found in provided reference. Seq: {seq}")
    return None


def MakeCoordinateLists(primerfile, referencefile, err_rate=0.1):
    keyl = ("LEFT", "PLUS", "POSITIVE", "FORWARD")
    keyr = ("RIGHT", "MINUS", "NEGATIVE", "REVERSE")

    seqs = list(SeqIO.parse(primerfile, "fasta"))

    def make_row(r):
        if any(o in r.id for o in keyl):
            return dict(Name=r.id, Seq=str(r.seq).upper(), orient="+")
        elif any(o in r.id for o in keyr):
            return dict(Name=r.id, Seq=str(r.seq).upper(), orient="-")
        print(
            f"Orientation for primer '{r.id}' not found. It is not included in the bed file."
        )
        return None

    df = pd.DataFrame(row for s in seqs if (row := make_row(s)) is not None)

    df["Seq_revcomp"] = df.apply(lambda x: Seq.reverse_complement(x["Seq"]), axis=1)

    ref_file = SeqIO.read(referencefile, "fasta")
    ref_seq = str(ref_file.seq)

    df["coords_fw"] = df.apply(
        lambda row: get_coords(row["Seq"], ref_seq, err_rate), axis=1
    )
    df["coords_rev"] = df.apply(
        lambda row: get_coords(row["Seq_revcomp"], ref_seq, err_rate), axis=1
    )
    df["ref"] = ref_file.name
    return df


def CoordinateListsToBed(df, outfile):
    with open(outfile, "w") as f:
        for i, row in df.iterrows():
            for start, stop in row.coords_fw:
                f.write(f"{row.ref}\t{start}\t{stop}\t{row.Name}\t.\t{row.orient}\n")


if __name__ == "__main__":
    import argparse

    args = argparse.ArgumentParser()

    args.add_argument(
        "--primers",
        metavar="File",
        type=str,
        help="Fasta file containing primers. Fasta headers should contain information about the orientation.",
        required=True,
    )

    args.add_argument(
        "--reference",
        metavar="File",
        type=str,
        help="Fasta file with the reference",
        required=True,
    )

    args.add_argument(
        "--output",
        metavar="File",
        type=str,
        help="The output BED file with coordinates of the primers.",
        required=True,
    )
    args.add_argument(
        "--primer-mismatch-rate",
        metavar="File",
        type=str,
        help="The fraction of mismatches a primer can have with respect to the reference.",
        default=0.1,
    )

    flags = args.parse_args()

    df = MakeCoordinateLists(
        flags.primers, flags.reference, flags["primer-mismatch-rate"]
    )

    CoordinateListsToBed(df, flags.output)
