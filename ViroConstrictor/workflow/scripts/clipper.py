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
    default=0,
)

arg.add_argument(
    "--max-aligned-length",
    metavar="FLOAT",
    type=float,
    required=False,
    help="Filter reads based on the maximum length of the aligned part of the read",
    default=0,
)

arg.add_argument(
    "--only-include-region",
    metavar="STRING",
    type=str,
    required=False,
    help="Only include reads where the aligned section of the read start and end within the specified region",
    default=None,
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
    """
    Split a CIGAR string into a list of tuples.

    Parameters:
    -----------
    cigar : str
        The CIGAR string to be split.

    Returns:
    --------
    list
        A list of tuples, where each tuple contains two elements:
        - The number of consecutive characters.
        - The character itself.

    Example:
    --------
    >>> split_cigar("3M1I2D")
    [(3, 'M'), (1, 'I'), (2, 'D')]
    """
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
    """
    Check if the given list of cigar tuples contains any spliced regions.

    Parameters
    ----------
    cigar_tuples : list
        A list of cigar tuples representing the alignment.

    Returns
    -------
    bool
        True if there are spliced regions, False otherwise.
    """
    return any(cigar_tuple[1] == "N" for cigar_tuple in cigar_tuples)


def get_largest_spliced_len(cigar_tuples: list) -> int:
    """
    Calculate the largest spliced length from a list of cigar tuples.

    Parameters:
    -----------
    cigar_tuples : list
        A list of cigar tuples representing the cigar string.

    Returns:
    --------
    int
        The largest spliced length found in the cigar tuples.

    Examples:
    ---------
    >>> cigar_tuples = [(10, 'M'), (5, 'N'), (20, 'M'), (15, 'N'), (30, 'M')]
    >>> get_largest_spliced_len(cigar_tuples)
    15
    """
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
    minimal_read_length = int(flags.min_aligned_length)
    maximum_read_length = int(
        reflength if flags.max_aligned_length == 0 else flags.max_aligned_length
    )

    include_region_start = (
        flags.only_include_region.split(":")[0]
        if flags.only_include_region is not None
        else None
    )
    include_region_end = (
        flags.only_include_region.split(":")[1]
        if flags.only_include_region is not None
        else None
    )

    for read in bamfile:
        read_start = read.query_alignment_start
        read_end = read.query_alignment_end
        ref_start = read.reference_start
        ref_end = read.reference_end
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
            if is_spliced(cigartuples) and (
                get_largest_spliced_len(cigartuples) > flags.spliced_length_threshold
            ):
                continue

        if read.query_alignment_length <= minimal_read_length:
            continue

        if read.query_alignment_length >= maximum_read_length:
            continue

        if include_region_start is not None and include_region_end is not None:
            if ref_start < int(include_region_start) or ref_end > int(
                include_region_end
            ):
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
