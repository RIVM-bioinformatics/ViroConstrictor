import sys

import pandas as pd
from Bio import SeqIO

ref, samplename, coverages, output = sys.argv[1:]

df = pd.read_csv(coverages, sep="\t", index_col=0)

for record in SeqIO.parse(ref, "fasta"):
    a = len(record.seq)

with open(output, "w") as f:
    f.write(
        f"{samplename}\t{100-(int(df[df < 1 ].count())/a)*100}\t{100-(int(df[df < 5 ].count())/a)*100}\t{100-(int(df[df < 10 ].count())/a)*100}\t{100-(int(df[df < 50 ].count())/a)*100}\t{100-(int(df[df < 100 ].count())/a)*100}\n"
    )
