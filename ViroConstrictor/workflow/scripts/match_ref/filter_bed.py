import os
import sys

import pandas as pd

bedfile, refdata, outfile, updatedstats = sys.argv[1:]

df = pd.read_csv(refdata, keep_default_na=False)

beddf = pd.read_csv(
    bedfile,
    sep="\t",
    names=["ref", "start", "stop", "name", "score", "strand"],
    index_col=False,
)

beddf = beddf.loc[beddf["ref"].isin(df["seqrecord_name"].tolist())].reset_index(
    drop=True
)

# replace the value of the "ref" column in beddf with the value of the "Reference" column in df, if the 'ref' value in beddf corresponds to the value of the 'seqrecord_name' column in df
beddf["ref"] = beddf["ref"].map(df.set_index("seqrecord_name")["Reference"])

beddf.to_csv(outfile, sep="\t", index=False, header=False)

df["Primer_file"] = os.path.abspath(outfile)
df.to_csv(updatedstats, index=False)
