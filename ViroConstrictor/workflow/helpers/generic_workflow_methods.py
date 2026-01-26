import AminoExtract
import numpy as np
import pandas as pd
from Bio import SeqIO, SeqRecord

from ViroConstrictor.workflow.helpers.directories import *
from ViroConstrictor.workflow.helpers.presets import get_preset_parameter


def get_reference_header(reffile):
    """
    #TODO: Docstring should be added to this function.
    #TODO: The docstrings in this file should also note where the function is used as they are explicit workflow helper functions.
    """
    return [record.id for record in SeqIO.parse(reffile, "fasta")]


def get_aminoacid_features(df):
    """
    #TODO: Docstring should be added to this function.
    #TODO: The docstrings in this file should also note where the function is used as they are explicit workflow helper functions.
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
    #TODO: Docstring should be added to this function.
    #TODO: The docstrings in this file should also note where the function is used as they are explicit workflow helper functions.
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
    #TODO: this docstring should be cleaned up and made more informative.
    #TODO: The docstrings in this file should also note where the function is used as they are explicit workflow helper functions.
    Read a FASTA file and return a list of SeqRecord objects.

    Parameters
    ----------
    fasta_file : str
        The path to the FASTA file to be read.

    Returns
    -------
    list[SeqRecord.SeqRecord]
        A list of SeqRecord objects representing the sequences in the FASTA file.

    """
    return list(SeqIO.parse(fasta_file, "fasta"))


def segmented_ref_groups(df: pd.DataFrame) -> pd.DataFrame:
    """
    #TODO: this docstring should be cleaned up and made more informative.
    #TODO: The docstrings in this file should also note where the function is used as they are explicit workflow helper functions.
    Filter required ref groups from main reference file.

    Parameters
    ----------
    df : pd.DataFrame
        A pandas DataFrame containing the reference file information.

    Returns
    -------
    pd.DataFrame
        A pandas DataFrame containing the filtered reference file information.

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


def get_features_all_samples(samples_df: pd.DataFrame) -> list[str]:
    """
    Extracts and returns a list of unique amino acid feature names across all samples  
    in the provided DataFrame.  

    Parameters  
    ----------  
    samples_df : pd.DataFrame  
        A pandas DataFrame containing data for multiple samples. Each row may contain  
        a column "AA_FEAT_NAMES" with a list or tuple of amino acid feature names.  

    Returns  
    -------  
    list of str  
        A list of unique amino acid feature names found across all samples in  
        ``samples_df``.
    """
    all_features = []
    for _, row in samples_df.iterrows():
        aa_feat_names = row.get("AA_FEAT_NAMES")
        if pd.notna(aa_feat_names) and isinstance(aa_feat_names, (list, tuple)):
            all_features.extend(aa_feat_names)
    return list(set(all_features))


def get_features_per_sample(sample: str, samples_df: pd.DataFrame) -> list[str]:
    """
    Extracts and returns a list of unique amino acid feature names for a sample.

    Parameters
    ----------
    sample : str
        Sample identifier.
    samples_df : pd.DataFrame
        A pandas DataFrame containing sample data.

    Returns
    -------
    list of str
        A list of unique amino acid feature names for the sample.
    """
    sample_rows = samples_df[samples_df["sample"] == sample]
    if sample_rows.empty:
        return []
    all_features = []
    for _, row in sample_rows.iterrows():
        aa_feat_names = row.get("AA_FEAT_NAMES")
        if pd.notna(aa_feat_names) and isinstance(aa_feat_names, (list, tuple)):
            all_features.extend(aa_feat_names)
    return list(set(all_features))

# Helper function to get all unique features for a virus
def get_features_per_virus(virus: str, samples_df: pd.DataFrame) -> list[str]:
    """
    Retrieve a list of unique amino acid feature names for a given virus from a DataFrame.

    Parameters
    ----------
    virus : str
        The name or identifier of the virus to filter the DataFrame.
    samples_df : pd.DataFrame
        A pandas DataFrame containing at least the columns "Virus" and "AA_FEAT_NAMES".
        The "AA_FEAT_NAMES" column should contain lists or tuples of feature names.

    Returns
    -------
    list of str
        A list of unique amino acid feature names associated with the specified virus.
        Returns an empty list if the virus is not found or no features are present.
    """
    virus_rows = samples_df[samples_df["Virus"] == virus]
    if virus_rows.empty:
        return []
    all_features = []
    for _, row in virus_rows.iterrows():
        aa_feat_names = row.get("AA_FEAT_NAMES")
        if pd.notna(aa_feat_names) and isinstance(aa_feat_names, (list, tuple)):
            all_features.extend(aa_feat_names)
    return list(set(all_features))