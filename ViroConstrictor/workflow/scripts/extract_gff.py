import sys

import gffpandas.gffpandas as gffpd

_, input, output, refID = sys.argv

gff = gffpd.read_gff3(input)

header = gff.header

gff.df = gff.df[gff.df.seq_id == refID]

gff.to_gff3(output)
