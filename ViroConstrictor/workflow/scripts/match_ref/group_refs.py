import os
import sys

import pandas as pd
from Bio import SeqIO
from rich import print

input_refs, input_stats, output_ref, output_stats, sample = sys.argv[1:]

input_refs = input_refs.split()
input_stats = input_stats.split()

pd.set_option("display.max_rows", None)
pd.set_option("display.max_columns", None)

seqrecords = []
for file in input_refs:
    seqrecords.extend(iter(SeqIO.parse(file, "fasta")))

df = pd.DataFrame(None)
# loop over the input stats files and append to the dataframe
for file in input_stats:
    tempdf = pd.read_csv(file, index_col=0)
    df = pd.concat([df, tempdf], ignore_index=True)

df["sample"] = sample

#TODO: refactor this to a more future-proof method. This is a temporary fix for the current issue and should be replaced with a method that doesn't have short-cut assumptions.
renamed_seqrecords = []
if len(seqrecords) > 1:
    # more than 1 secrecord assumes segmented mode.
    # We will therefore assume the reference headers are formatted to the segment-mode format on docs-site.
    # we also now have to rename the records for proper future processing of all segments
    # rename the seqrecords where record.id becomes the first item in the record.description (split on "|")
    for record in seqrecords:
        record.id = record.description.split()[1].split("|")[0]
        record.description = " ".join(
            [record.name, " ".join(record.description.split(" ")[1:])]
        )
        renamed_seqrecords.append(record)
else:
    # only 1 seqrecord assumes non-segmented mode
    # we will however then assume the reference is formatted like an NCBI nuccore accession
    # > record.id = nuccore accession || record.description = pathogen description
    renamed_seqrecords.append(seqrecords[0])

# add the seqrecord information to the dataframe
for record in renamed_seqrecords:
    df.loc[df["Reference"].str.contains(record.id), "seqrecord_id"] = record.id
    df.loc[df["Reference"].str.contains(record.id), "seqrecord_description"] = (
        record.description
    )
    df.loc[df["Reference"].str.contains(record.id), "seqrecord_name"] = record.name
    df.loc[df["Reference"].str.contains(record.id), "seqrecord_seq"] = str(record.seq)

# replace the "Reference" column with the "seqrecord_id" column
df["Reference"] = df["seqrecord_id"]
df["Reference_file"] = os.path.abspath(output_ref)

# write the dataframe to a csv
df.to_csv(output_stats, index=False)

SeqIO.write(renamed_seqrecords, output_ref, "fasta")
