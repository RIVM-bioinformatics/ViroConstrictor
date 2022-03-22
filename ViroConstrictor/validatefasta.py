# pylint: disable=C0103

"""
Basic functions to see if a fasta is valid
"""

import re

from Bio import SeqIO


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
        return not any(ContainsAmbiguities(str(record.seq)) for record in SeqIO.parse(inputfile, "fasta"))
    return False

def IsValidFasta(inputfile):
    if inputfile == "NONE":
        return True
    results = [
        ContainsSpecials(str(record.seq))
        for record in SeqIO.parse(inputfile, "fasta")
    ]

    if any(results):
        return False
    return True
