import sys

import AminoExtract

_, input, output, refID = sys.argv

gff = AminoExtract.read_gff(input)

header = gff.header
gff.df = gff.df[gff.df.seqid == refID]

with open(output, "w") as f:
    f.write(header)
    f.write(gff.df.to_csv(sep="\t", index=False, header=None))
