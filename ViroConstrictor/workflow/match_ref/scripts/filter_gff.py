import os
import sys

import AminoExtract
import pandas as pd

refdata, gfffile, outfile, updatedstats = sys.argv[1:]

gff = AminoExtract.read_gff(gfffile, split_attributes=False)

df = pd.read_csv(refdata, keep_default_na=False)

seqrecordNames = df["seqrecord_name"].tolist()
gff.df = gff.df.loc[gff.df["seqid"].isin(seqrecordNames)].reset_index(drop=True)

df["Feat_file"] = os.path.abspath(outfile)

# replace the value of the "seqid" column in gffdf with the value of the "Reference" column in df, if the 'seqid' value in gffdf corresponds to the value of the 'seqrecord_name' column in df
gff.df["seqid"] = gff.df["seqid"].map(df.set_index("seqrecord_name")["Reference"])

feats_to_write = gff.df.to_csv(sep="\t", index=False, header=False)


with open(outfile, "w") as f:
    f.write(gff.header)
    f.write(feats_to_write)

df.to_csv(updatedstats, index=False)
