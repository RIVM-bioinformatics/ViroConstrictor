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


def IsValidFasta(inputfile):
    if inputfile == "NONE":
        return True
    results = []
    for record in SeqIO.parse(inputfile, "fasta"):
        results.append(ContainsSpecials(str(record.seq)))

    if any(results) is True:
        return False
    return True
