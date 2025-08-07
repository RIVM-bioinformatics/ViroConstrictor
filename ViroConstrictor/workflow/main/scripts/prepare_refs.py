import sys

from Bio import SeqIO

_, input, output, reference_id = sys.argv

for record in SeqIO.parse(str(input), "fasta"):
    if reference_id in record.id:
        record.seq = record.seq.upper()
        SeqIO.write(record, str(output), "fasta")
