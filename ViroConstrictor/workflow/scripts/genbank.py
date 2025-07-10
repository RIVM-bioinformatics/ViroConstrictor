import sys
from pathlib import Path

from BCBio import GFF
from Bio import SeqIO


def open_genbank_file(file_path: Path) -> list[SeqIO.SeqRecord]:
    """
    Opens a GenBank file and returns a list of SeqIO.SeqRecord objects.

    Args:
        file_path (Path): The path to the GenBank file.

    Returns:
        List of SeqIO.SeqRecord objects.
    """
    try:
        return list(SeqIO.parse(file_path, "genbank"))
    except Exception as e:
        raise ValueError(f"Error opening GenBank file: {e}")


def write_target(file_path: Path, records: list[SeqIO.SeqRecord]) -> None:
    organisms: list[str | None] = [
        record.annotations.get("organism", None) for record in records
    ]
    organisms = [
        org.split("(")[0].strip().replace(" ", "_") for org in organisms if org
    ]
    if not all(org == organisms[0] for org in organisms):
        raise ValueError(
            "Not all GenBank records have the same organism annotation.\n"
            "Either edit the GenBank file to have the same organism annotation for all records (strains don't count),\n"
            "or manually provide a target organism name using the --target option."
        )
    org = organisms[0]

    target_file_path = file_path.parent / f"{org}"

    with open(target_file_path, "w", encoding="utf-8") as target_file:
        # write all records as target to a single file
        target_file.write(f"{org}\n")


def split_genbank(file_path: Path, include_target: bool = False) -> None:
    """Splits a GenBank file into a reference fasta, a features file and possibly a target file."""
    records = open_genbank_file(file_path)
    with open(file_path.with_suffix(".fasta"), "w", encoding="utf-8") as fasta_file:
        # write all records as fasta to a single file
        for record in records:
            SeqIO.write(record, fasta_file, "fasta")

            # SeqIO.write(record, gff_file, "gff")
    with open(file_path.with_suffix(".gff"), "w", encoding="utf-8") as gff_file:
        GFF.write(records, gff_file)

    if include_target:
        write_target(file_path, records)


def main(args: list[str] | None = None) -> None:
    """Main function to run the script."""
    if args is None:
        args = sys.argv[1:]

    if len(args) < 1:
        print("Usage: genbank.py <genbank_file> [--target]")
        sys.exit(1)

    file_path = Path(args[0])
    include_target = "--target" in args

    split_genbank(file_path, include_target)


if __name__ == "__main__":
    main()
