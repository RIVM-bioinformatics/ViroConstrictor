import sys

import pandas as pd
from Bio import SeqIO

inputcount, inputref, filtref, filtcount = sys.argv[1:]

# Read in dataframe
df = pd.read_csv(inputcount, index_col=0)

# sort by "Mapped Reads" column
df = df.sort_values(by="Mapped Reads", ascending=False)
# keep only the first row
df = df.iloc[[0]]

# get the reference name
refname = df["Reference"].values[0]

# read the reference fasta file and filter for the reference name
for record in SeqIO.parse(inputref, "fasta"):
    if record.id == refname:
        SeqIO.write(record, filtref, "fasta")

df.to_csv(filtcount)
