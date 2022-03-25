from Bio import SeqIO
import sys

_, input_ref, output_file, refID_wildcard = sys.argv


if __name__ == "__main__":
    for record in SeqIO.parse(input_ref, "fasta"):
        if refID_wildcard.split('~')[1] in record.id:
            SeqIO.write(record, output_file, "fasta")