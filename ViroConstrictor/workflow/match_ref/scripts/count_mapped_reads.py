import sys

import pandas as pd
import pysam

inputfile, outputfile = sys.argv[1:]

# Open BAM file
bamfile = pysam.AlignmentFile(inputfile, "rb")

# Loop over reference headers and count mapped reads, mismatches, and sequence identity
data = []
for ref in bamfile.references:
    mapped_reads = 0
    total_mismatches = 0
    total_identity = 0
    for read in bamfile.fetch(ref):
        if not read.is_unmapped:
            mapped_reads += 1
            total_mismatches += read.get_tag("NM")
            total_identity += 1 - (read.get_tag("NM") / read.query_alignment_length)
    if mapped_reads > 0:
        avg_mismatches = total_mismatches / mapped_reads
        avg_identity = total_identity / mapped_reads
    else:
        avg_mismatches = 0
        avg_identity = 0
    data.append(
        {
            "Reference": ref,
            "Mapped Reads": mapped_reads,
            "Avg. Mismatches per Read": avg_mismatches,
            "Avg. Sequence Identity": avg_identity,
        }
    )

# Create pandas dataframe
df = pd.DataFrame(data)

# output to feather file
df.to_csv(outputfile)
