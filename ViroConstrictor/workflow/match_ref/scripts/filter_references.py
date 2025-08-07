import sys

from Bio import SeqIO

_, input, output, wildcard_segment = sys.argv

records_to_keep = []
for record in SeqIO.parse(str(input), "fasta"):
    if wildcard_segment != "None":  # wildcard_segment is a string instead of a NoneType
        if record.description.split(" ")[1].split("|")[0] == wildcard_segment:
            records_to_keep.append(record)
    else:
        records_to_keep.append(record)
SeqIO.write(records_to_keep, str(output), "fasta")
