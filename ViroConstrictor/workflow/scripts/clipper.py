import argparse

import pysam

arg = argparse.ArgumentParser()

arg.add_argument(
    "--input",
    metavar="File",
    help="BAM file of the first minimap2 alignment run, sorted and indexed.",
    type=str,
    required=True,
)

arg.add_argument(
    "--output",
    metavar="File",
    help="File with cleaned fastq reads",
    type=argparse.FileType("w"),
    required=True,
)

arg.add_argument(
    "--exclude-spliced",
    action="store_true",
    default=False,
    help="Exclude spliced reads from the output",
    required=False,
)

arg.add_argument(
    "--spliced-length-threshold",
    metavar="INT",
    type=int,
    default=0,
    help="When excluding spliced reads, only exclude reads with a spliced length above this threshold",
    required=False,
)

arg.add_argument(
    "--min-aligned-length",
    metavar="FLOAT",
    type=float,
    required=False,
    help="Filter reads based on the minimum length of the aligned part of the read",
    default=1,
)

arg.add_argument(
    "--threads",
    metavar="Number",
    help="Number of threads that can be used for decompressing/compressing the BAM file",
    default=1,
    type=int,
    required=False,
)

flags = arg.parse_args()


def split_cigar(cigar: str) -> list:
    cigar_tuples = []
    current_number = ""
    for char in cigar:
        if char.isdigit():
            current_number += char
        else:
            cigar_tuples.append((int(current_number), char))
            current_number = ""
    return cigar_tuples


def is_spliced(cigar_tuples: list) -> bool:
    return any(cigar_tuple[1] == "N" for cigar_tuple in cigar_tuples)


def get_largest_spliced_len(cigar_tuples: list) -> int:
    largest_spliced_len = 0
    current_spliced_len = 0
    for cigar_tuple in cigar_tuples:
        if cigar_tuple[1] == "N":
            current_spliced_len += cigar_tuple[0]
        else:
            if current_spliced_len > largest_spliced_len:
                largest_spliced_len = current_spliced_len
            current_spliced_len = 0
    return largest_spliced_len


with flags.output as fileout:
    bamfile = pysam.AlignmentFile(flags.input, "rb", threads=flags.threads)

    reflength = bamfile.lengths[0]
    minimal_read_length = int(reflength * flags.min_aligned_length)

    for read in bamfile:

        read_start = read.query_alignment_start
        read_end = read.query_alignment_end

        trimmed_seq = read.query_alignment_sequence
        trimmed_qual = read.qual[read_start:read_end]

        if read.is_reverse == True:
            complement = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
            bases = list(trimmed_seq)
            bases = [complement[base] for base in bases]
            trimmed_seq = "".join(bases)
            trimmed_seq = trimmed_seq[::-1]
            trimmed_qual = trimmed_qual[::-1]

        if len(trimmed_seq) == 0:
            continue

        if flags.exclude_spliced:
            cigartuples = split_cigar(read.cigarstring)
            if is_spliced(cigartuples):
                if (
                    get_largest_spliced_len(cigartuples)
                    < flags.spliced_length_threshold
                ):
                    continue

        if read.query_alignment_length <= minimal_read_length:
            continue

        fileout.write(
            "@"
            + str(read.query_name)
            + "\n"
            + str(trimmed_seq)
            + "\n"
            + "+"
            + "\n"
            + str(trimmed_qual)
            + "\n"
        )

fileout.close()
