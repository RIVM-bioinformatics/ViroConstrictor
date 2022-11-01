import pandas as pd
import sys

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

# change the INFO column to only contain the depth-numbers, remove the keys
df["INFO"] = df["INFO"].str.split("=", expand=True)[1].str.split(";", expand=True)[0]

# drop the unnecessary columns in the dataframe
df.drop(["ID", "QUAL", "FILTER"], axis=1, inplace=True)

# Rename the columns to match the format
df.rename(columns={
    "CHROM": "Reference_ID",
    "POS": "Position",
    "REF": "Reference_Base",
    "ALT": "Variant_Base",
    "INFO": "Depth"
    }, inplace=True)

# insert the samplename as a column
df.insert(loc=0, column="Sample", value=samplename)

# write the file as a TSV to disk
df.to_csv(output_tsv, sep="\t", index=False, header=False)