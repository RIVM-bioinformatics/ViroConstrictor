import inspect
import re

import AminoExtract
import numpy as np
import pandas as pd
from Bio import SeqIO, SeqRecord

from ViroConstrictor.workflow.helpers.directories import *
from ViroConstrictor.workflow.helpers.presets import get_preset_parameter


def get_reference_header(reffile):
    """
    Extract reference sequence identifiers from a FASTA file.

    Used in main/workflow.smk to retrieve reference IDs for wildcard expansion.

    Parameters
    ----------
    reffile : str or Path
        Path to the FASTA reference file.

    Returns
    -------
    list[str]
        List of sequence identifiers (record IDs) from the FASTA file.
    """
    return [record.id for record in SeqIO.parse(reffile, "fasta")]


def get_aminoacid_features(df):
    """
    Extract amino acid feature names from reference GFF files for each sample.

    Used in main/workflow.smk to populate the AA_FEAT_NAMES column based on preset
    parameters and GFF feature types. Returns NaN for samples without features.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with columns FEATURES (path to GFF file or "NONE"), REFERENCE (path
        to reference sequence), PRESET (preset name), and RefID (reference identifier).

    Returns
    -------
    pd.DataFrame
        DataFrame with added/updated AA_FEAT_NAMES column containing tuples of feature
        names or NaN where features are unavailable.
    """
    records = df.to_dict(orient="records")

    for rec in records:
        if rec["FEATURES"] != "NONE":
            AA_dict = AminoExtract.get_feature_name_attribute(
                input_gff=str(rec["FEATURES"]),
                input_seq=str(rec["REFERENCE"]),
                feature_types=get_preset_parameter(rec["PRESET"], "AminoExtract_FeatureType", "GLOBAL"),
            )
            if AA_dict:
                for k, v in AA_dict.items():
                    if k == rec["RefID"]:
                        rec["AA_FEAT_NAMES"] = tuple(v)
            else:
                rec["AA_FEAT_NAMES"] = np.nan
        else:
            rec["AA_FEAT_NAMES"] = np.nan
    return pd.DataFrame.from_records(records)


def list_aminoacid_result_outputs(samples_df):
    """
    Generate list of amino acid output file paths for samples with AA features.

    Used in main/workflow.smk to construct expected output rules for amino acid
    sequences. Generates filenames based on Virus, RefID, and feature names.

    Parameters
    ----------
    samples_df : pd.DataFrame
        DataFrame with columns Virus (virus name), RefID (reference ID), and
        AA_FEAT_NAMES (tuple or NaN of feature names).

    Returns
    -------
    list[str]
        List of unique output file paths in format ``{res}Virus~{Virus}/RefID~{RefID}/{amino}{aa}.faa``,
        or empty list if no samples have AA features.
    """
    aminoacid_features = []
    for x in samples_df.to_dict(orient="records"):
        Virus = x["Virus"]
        RefID = x["RefID"]
        features = x["AA_FEAT_NAMES"]
        if not isinstance(features, float):
            aminoacid_features.extend([f"{res}Virus~{Virus}/RefID~{RefID}/{amino}{aa}.faa" for aa in features])
    return list(set(aminoacid_features))


def read_fasta(fasta_file: str) -> list[SeqRecord.SeqRecord]:
    """
    Read sequences from a FASTA file.

    Used in main/workflow.smk and match_ref/workflow.smk for parsing reference sequences.

    Parameters
    ----------
    fasta_file : str
        Path to the FASTA file to be read.

    Returns
    -------
    list[SeqRecord.SeqRecord]
        List of SeqRecord objects representing the sequences in the FASTA file.
    """
    return list(SeqIO.parse(fasta_file, "fasta"))


def segmented_ref_groups(df: pd.DataFrame) -> pd.DataFrame:
    """
    Filter and segment reference groups from a reference metadata DataFrame.

    Used in match_ref/workflow.smk to prepare reference information for segmented
    viruses. Extracts segment identifiers from FASTA sequence descriptions and
    removes references with fewer than 2 segments.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing reference metadata with columns SEGMENTED (bool),
        REFERENCE (path to FASTA), and segment (to be populated).

    Returns
    -------
    pd.DataFrame
        Updated DataFrame with segment column populated as sets of segment names
        for segmented viruses, or {"None"} for non-segmented references.
        Rows with single segment are removed.

    Notes
    -----
    Segment identifiers are extracted from FASTA description field (second word,
    before first "|"). Sets are flattened after extraction.
    """
    for index, row in df.iterrows():
        # if the value in the "SEGMENTED" column is False then place a None string in the segment column.
        # Ensure that the value is a string and not a NoneType.
        if not row["SEGMENTED"]:
            df.at[index, "segment"] = {"None"}
            continue
        refs = read_fasta(row["REFERENCE"])
        unique_groups = {x.description.split(" ")[1].split("|")[0] for x in refs}
        if len(unique_groups) < 2:
            df.drop(index, inplace=True)
            continue
        df.at[index, "segment"] = [unique_groups]
    # ensure the segment column is a flat set of strings and not a list with a set inside.
    df["segment"] = df["segment"].apply(lambda x: x.pop() if isinstance(x, list) else x)
    # TODO: Check this later to see if this can be done in the first pass (when the values are being set in the first place) instead of fixing afterwards.
    # TODO: Not entirely sure why the values are being generated as a list with a set inside one time and just as a normal set the other time.
    return df


# TODO: all functions below this lack proper unit tests. This is acceptable (at time of writing) for the current release. However, proper unit tests should be added in the near future.
def get_features_all_samples(samples_df: pd.DataFrame) -> list[str]:
    """
    Extract unique amino acid feature names across all samples.

    Used in main/workflow.smk to generate aggregate amino acid output targets.

    Parameters
    ----------
    samples_df : pd.DataFrame
        DataFrame containing sample data with AA_FEAT_NAMES column (list, tuple,
        or NaN of amino acid feature names per sample).

    Returns
    -------
    list[str]
        Sorted list of unique amino acid feature names found across all samples.
    """
    all_features = []
    for _, row in samples_df.iterrows():
        aa_feat_names = row.get("AA_FEAT_NAMES")
        if isinstance(aa_feat_names, (list, tuple)):
            all_features.extend(aa_feat_names)
    return sorted(list(set(all_features)))


def get_features_per_sample(sample: str, samples_df: pd.DataFrame) -> list[str]:
    """
    Extract unique amino acid feature names for a specific sample.

    Used in main/workflow.smk to generate per-sample amino acid output targets.

    Parameters
    ----------
    sample : str
        Sample identifier to filter by.
    samples_df : pd.DataFrame
        DataFrame containing sample data with "sample" and "AA_FEAT_NAMES" columns.

    Returns
    -------
    list[str]
        Sorted list of unique amino acid feature names for the sample, or empty list
        if sample not found or has no features.
    """
    sample_rows = samples_df[samples_df["sample"] == sample]
    if sample_rows.empty:
        return []
    all_features = []
    for _, row in sample_rows.iterrows():
        aa_feat_names = row.get("AA_FEAT_NAMES")
        if isinstance(aa_feat_names, (list, tuple)):
            all_features.extend(aa_feat_names)
    return sorted(list(set(all_features)))


# Helper function to get all unique features for a virus
def get_features_per_virus(virus: str, samples_df: pd.DataFrame) -> list[str]:
    """
    Extract unique amino acid feature names for a specific virus.

    Used in main/workflow.smk to generate per-virus amino acid output targets.

    Parameters
    ----------
    virus : str
        Virus identifier to filter by.
    samples_df : pd.DataFrame
        DataFrame containing sample data with "Virus" and "AA_FEAT_NAMES" columns.
        AA_FEAT_NAMES should contain lists or tuples of feature names.

    Returns
    -------
    list[str]
        Sorted list of unique amino acid feature names associated with the virus,
        or empty list if virus not found or has no features.
    """
    virus_rows = samples_df[samples_df["Virus"] == virus]
    if virus_rows.empty:
        return []
    all_features = []
    for _, row in virus_rows.iterrows():
        aa_feat_names = row.get("AA_FEAT_NAMES")
        if isinstance(aa_feat_names, (list, tuple)):
            all_features.extend(aa_feat_names)
    return sorted(list(set(all_features)))

def get_rule_name() -> str:
    """
    Return the name of the closest Snakemake rule above the current line.

    This function inspects the Python call stack to access the workflow
    source code stored in the calling frame and determines which rule
    decorator appears closest above the current execution line.

    It searches for rule decorators of the form::

        @workflow.rule(name="rule_name", lineno=N)

    and returns the rule whose ``lineno`` is the largest value less than
    or equal to the current execution line.

    Returns
    -------
    str
        The name of the closest matching rule. If the rule name cannot
        be inferred, the fallback value ``"unknown_rule"`` is returned.

    Notes
    -----
    The function inspects two frames up the call stack to access the rule definition context. 
    The frame reference is deleted afterward to avoid reference cycles when using the `inspect` module.
    """
    fallback = "unknown_rule"
    frame = inspect.currentframe()
    try:
        caller = frame.f_back if frame else None
        rule_frame = caller.f_back if caller else None
        code = getattr(rule_frame.f_locals, "get", lambda k, d=None: d)("code", "") if rule_frame else ""
        lineno = getattr(caller, "f_lineno", None) if caller else None

        if not code or not isinstance(lineno, int):
            return fallback

        # Name of the closest matching rule
        return max(
            ((m.group(1), int(m.group(2))) for m in re.finditer(
                r"@workflow\.rule\(name=['\"]([^'\"]+)['\"],\s*lineno=(\d+)", code)
             if int(m.group(2)) <= lineno),
            default=(fallback, 0),
            key=lambda x: x[1]
        )[0]

    finally:
        del frame
