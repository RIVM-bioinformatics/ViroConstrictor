import os
import sys

import pandas as pd

refdata, gfffile, outfile, updatedstats = sys.argv[1:]

gffheader = ""
with open(gfffile, "r") as f:
    for line in f:
        if line.startswith("#"):
            gffheader += line
        else:
            break

gffdf = pd.read_csv(
    gfffile,
    sep="\t",
    comment="#",
    names=[
        "seqid",
        "source",
        "type",
        "start",
        "end",
        "score",
        "strand",
        "phase",
        "attributes",
    ],
)

df = pd.read_csv(refdata, keep_default_na=False)

seqrecordNames = df["seqrecord_name"].tolist()
gffdf = gffdf.loc[gffdf["seqid"].isin(seqrecordNames)].reset_index(drop=True)

df["Feat_file"] = os.path.abspath(outfile)

# replace the value of the "seqid" column in gffdf with the value of the "Reference" column in df, if the 'seqid' value in gffdf corresponds to the value of the 'seqrecord_name' column in df
gffdf["seqid"] = gffdf["seqid"].map(df.set_index("seqrecord_name")["Reference"])

feats_to_write = gffdf.to_csv(sep="\t", index=False, header=False)


with open(outfile, "w") as f:
    f.write(gffheader)
    f.write(feats_to_write)

df.to_csv(updatedstats, index=False)
