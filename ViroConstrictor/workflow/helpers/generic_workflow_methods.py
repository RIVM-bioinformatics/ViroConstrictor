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
