import sys

import pandas as pd

_, input_vcf, output_tsv, samplename = sys.argv

# read the file
df = pd.read_csv(
    input_vcf,
    sep="\t",
    comment="#",
    names=["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"],
)

# drop all rows where ALT is a non-determined base (N)
df = df[df.ALT != "N"]

if df.empty is True:
    df.to_csv(output_tsv, sep="\t", index=False, header=False)
    sys.exit(0)

# change the INFO column to only contain the depth-numbers, remove the keys
df["INFO"] = df["INFO"].str.split("=", expand=True)[1].str.split(";", expand=True)[0]

# drop the unnecessary columns in the dataframe
df.drop(["ID", "QUAL", "FILTER"], axis=1, inplace=True)

# insert the samplename as a column
df.insert(loc=0, column="Sample", value=samplename)

# write the file as a TSV to disk
df.to_csv(output_tsv, sep="\t", index=False, header=False)
