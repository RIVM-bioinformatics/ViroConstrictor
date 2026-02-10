import copy
from typing import Literal

import pandas as pd

from ViroConstrictor.logging import log
from ViroConstrictor.parser import CLIparser
from ViroConstrictor.runreport import WriteReport
from ViroConstrictor.scheduler import Scheduler
from ViroConstrictor.workflow_executor import run_snakemake_workflow


def replace_sets_to_singular_values(df: pd.DataFrame, columns: list) -> pd.DataFrame:
    """
    Replace sets in specified columns with their singular values.

    Parameters
    ----------
    df : pd.DataFrame
        The DataFrame containing the columns to be processed.
    columns : list
        The list of column names to be processed.

    Returns
    -------
    pd.DataFrame
        The DataFrame with the specified columns' sets replaced with their singular values.
    """
    for col in columns:
        df[col] = df[col].apply(lambda x: list(x)[0])
    return df


def replacement_merge_dataframe_on_cols(
    original_df: pd.DataFrame,
    override_df: pd.DataFrame,
    cols_left: list,
    cols_right: list,
) -> pd.DataFrame:
    """
    Merge two dataframes based on a common column and replace values in the original dataframe with values from the override dataframe.

    Parameters
    ----------
    original_df : pd.DataFrame
        The original dataframe to be merged and updated.
    override_df : pd.DataFrame
        The dataframe containing the updated values to be merged into the original dataframe.
    cols_left : list
        The list of column names in the original dataframe to be merged.
    cols_right : list
        The list of column names in the override dataframe to be merged.

    Returns
    -------
    pd.DataFrame
        The merged dataframe with updated values from the override dataframe.
    """

    # set sample columns to str type to avoid issues with merging
    original_df["SAMPLE"] = original_df["SAMPLE"].astype(str)
    override_df["sample"] = override_df["sample"].astype(str)

    for replacement_columns in zip(cols_left, cols_right):
        original_df[replacement_columns[0]] = original_df.apply(
            lambda x, replacement_columns=replacement_columns: (
                override_df[replacement_columns[1]][override_df["sample"] == x["SAMPLE"]].values[0]
                if x["SAMPLE"] in override_df["sample"].values and x[replacement_columns[0]] != "NONE"
                else x[replacement_columns[0]]
            ),
            axis=1,
        )
    return original_df


def process_match_ref(parsed_inputs: CLIparser, scheduler: Scheduler) -> CLIparser:
    """
    Process the match_ref step of the ViroConstrictor pipeline.

    This function runs the snakemake pipeline for the match_ref step using the parsed inputs. If the pipeline fails, it writes a report and exits with status code 1. If the pipeline succeeds, it reads the resulting dataframe from the match_ref_results.pkl file and extracts the columns "sample", "Reference", "Reference_file", "Primer_file" and "Feat_file". It then groups the dataframe by sample and aggregates the values in each column into a set. The resulting dataframe is then modified to replace all the sets with just one value to only that value (without the set) except for the Reference column. Finally, the values for the columns "REFERENCE", "PRIMERS" and "FEATURES" in parsed_inputs.samples_df are replaced with the values in columns "Reference_file", "Primer_file" and "Feat_file" in the modified dataframe where the sample names match.

    Parameters
    ----------
    parsed_inputs : CLIparser
        The parsed inputs for the ViroConstrictor pipeline.

    Returns
    -------
    CLIparser
        The modified parsed inputs with updated values for the "REFERENCE", "PRIMERS" and "FEATURES" columns.
    """
    inputs_obj_match_ref = copy.deepcopy(parsed_inputs)

    log.info(f"{'='*20} [bold orange_red1] Starting Match-reference process [/bold orange_red1] {'='*20}")

    status, used_workflow_config = run_snakemake_workflow(inputs_obj_match_ref, stage="MR", scheduler=scheduler)

    workflow_state: Literal["Failed", "Success"] = "Failed" if status is False else "Success"
    if status is False:
        WriteReport(
            inputs_obj_match_ref.workdir,
            inputs_obj_match_ref.input_path,
            inputs_obj_match_ref.exec_start_path,
            inputs_obj_match_ref.user_config,
            used_workflow_config.resource_settings,
            used_workflow_config.output_settings,
            inputs_obj_match_ref,
            workflow_state,
        )
        exit(1)

    if used_workflow_config.output_settings.dryrun is True:
        return parsed_inputs

    # status is True, this means that a file should exist at {workdir}/data/match_ref_results.pkl
    df = pd.read_pickle(f"{inputs_obj_match_ref.workdir}/data/match_ref_results.pkl")

    # ensure the "Primer_file" and "Feat_file" columns exist in df. If they don't, create them and set their values to "NONE"
    if "Primer_file" not in df.columns:
        df["Primer_file"] = "NONE"
    if "Feat_file" not in df.columns:
        df["Feat_file"] = "NONE"

    # keep only columns "sample", "Reference", "Reference_file", "Primer_file" and "Feat_file" from df
    filt_df = df[["sample", "Reference", "Reference_file", "Primer_file", "Feat_file"]].copy()

    # replace any cell in columns "Primer_file" and "Feat_file" that contains a NaN or None value with the string "NONE"
    filt_df.loc[:, "Primer_file"] = filt_df["Primer_file"].fillna("NONE")
    filt_df.loc[:, "Feat_file"] = filt_df["Feat_file"].fillna("NONE")

    imploded_df = filt_df.groupby("sample").agg(set).reset_index()

    # replace all the sets with just one value to only that value (without the set) except for the Reference column
    imploded_df = replace_sets_to_singular_values(imploded_df, ["Reference_file", "Primer_file", "Feat_file"])

    # replace the values for the columns "REFERENCE", "PRIMERS" and "FEATURES" in parsed_inputs.samples_df with the values in columns "Reference_file", "Primer_file" and "Feat_file" in imploded_df where the sample names match.
    # if the sample names don't match, the value in the column should be unchanged from the original value
    parsed_inputs.samples_df = replacement_merge_dataframe_on_cols(
        parsed_inputs.samples_df,
        imploded_df,
        ["REFERENCE", "PRIMERS", "FEATURES"],
        ["Reference_file", "Primer_file", "Feat_file"],
    )

    parsed_inputs.samples_df.set_index("SAMPLE", inplace=True)
    parsed_inputs.samples_dict = parsed_inputs.samples_df.to_dict(orient="index")

    return parsed_inputs
