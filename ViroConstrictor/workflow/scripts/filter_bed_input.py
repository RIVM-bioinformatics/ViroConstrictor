import sys

import pandas as pd

_, inputfile, outputfile, refID = sys.argv

df = pd.read_csv(
    inputfile, sep="\t", names=["ref", "start", "stop", "name", "score", "strand"], index_col=False
)

df = df[df.ref == refID]

df.to_csv(outputfile, sep="\t", index=False, header=False)
