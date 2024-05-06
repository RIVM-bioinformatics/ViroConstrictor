import copy
import sys
from typing import Literal

import pandas as pd
from snakemake import snakemake

import ViroConstrictor.logging
from ViroConstrictor.logging import log
from ViroConstrictor.parser import CLIparser
from ViroConstrictor.runconfigs import GetSnakemakeRunDetails, WriteYaml
from ViroConstrictor.runreport import WriteReport
from ViroConstrictor.workflow.containers import (
    construct_container_bind_args,
    download_containers,
)


def run_snakemake(
    inputs_obj: CLIparser, snakemakedetails: GetSnakemakeRunDetails
) -> bool:
    """
    Run Snakemake pipeline.

    Parameters
    ----------
    inputs_obj : CLIparser
        An instance of the CLIparser class containing the parsed command-line arguments.
    snakemakedetails : GetSnakemakeRunDetails
        An instance of the GetSnakemakeRunDetails class containing the details of the Snakemake run.

    Returns
    -------
    bool
        True if the Snakemake pipeline ran successfully, False otherwise.
    """
    if inputs_obj.user_config["COMPUTING"]["compmode"] == "local":
        return snakemake(
            snakefile=inputs_obj.match_ref_snakefile,
            workdir=inputs_obj.workdir,
            cores=snakemakedetails.snakemake_run_conf["cores"],
            use_conda=snakemakedetails.snakemake_run_conf["use-conda"],
            conda_frontend="mamba",
            use_singularity=snakemakedetails.snakemake_run_conf["use-singularity"],
            singularity_args=construct_container_bind_args(inputs_obj.samples_dict),
            jobname=snakemakedetails.snakemake_run_conf["jobname"],
            latency_wait=snakemakedetails.snakemake_run_conf["latency-wait"],
            dryrun=snakemakedetails.snakemake_run_conf["dryrun"],
            configfiles=[
                WriteYaml(
                    snakemakedetails.snakemake_run_parameters,
                    f"{inputs_obj.workdir}/config/run_params_MR.yaml",
                ),
                WriteYaml(
                    snakemakedetails.snakemake_run_conf,
                    f"{inputs_obj.workdir}/config/run_configs_MR.yaml",
                ),
            ],
            restart_times=snakemakedetails.snakemake_run_conf["restart-times"],
            keepgoing=snakemakedetails.snakemake_run_conf["keep-going"],
            quiet=["all"],  # type: ignore
            log_handler=[
                ViroConstrictor.logging.snakemake_logger(logfile=inputs_obj.logfile),
            ],
            printshellcmds=snakemakedetails.snakemake_run_conf["printshellcmds"],
            scheduler=snakemakedetails.snakemake_run_conf["scheduler"],
        )

    return snakemake(
        snakefile=inputs_obj.match_ref_snakefile,
        workdir=inputs_obj.workdir,
        cores=snakemakedetails.snakemake_run_conf["cores"],
        use_conda=snakemakedetails.snakemake_run_conf["use-conda"],
        conda_frontend="mamba",
        use_singularity=snakemakedetails.snakemake_run_conf["use-singularity"],
        singularity_args=construct_container_bind_args(inputs_obj.samples_dict),
        jobname=snakemakedetails.snakemake_run_conf["jobname"],
        latency_wait=snakemakedetails.snakemake_run_conf["latency-wait"],
        drmaa=snakemakedetails.snakemake_run_conf["drmaa"],
        drmaa_log_dir=snakemakedetails.snakemake_run_conf["drmaa-log-dir"],
        dryrun=snakemakedetails.snakemake_run_conf["dryrun"],
        configfiles=[
            WriteYaml(
                snakemakedetails.snakemake_run_parameters,
                f"{inputs_obj.workdir}/config/run_params_MR.yaml",
            ),
            WriteYaml(
                snakemakedetails.snakemake_run_conf,
                f"{inputs_obj.workdir}/config/run_configs_MR.yaml",
            ),
        ],
        restart_times=snakemakedetails.snakemake_run_conf["restart-times"],
        keepgoing=snakemakedetails.snakemake_run_conf["keep-going"],
        quiet=["all"],  # type: ignore
        log_handler=[
            ViroConstrictor.logging.snakemake_logger(logfile=inputs_obj.logfile),
        ],
        printshellcmds=snakemakedetails.snakemake_run_conf["printshellcmds"],
        scheduler=snakemakedetails.snakemake_run_conf["scheduler"],
    )


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


def get_snakemake_run_details(inputs_obj: CLIparser) -> GetSnakemakeRunDetails:
    """
    Filter the samples dataframe to only include samples that require matching to a reference sequence.
    Set the index of the filtered dataframe to the "SAMPLE" column.
    Convert the filtered dataframe to a dictionary with "SAMPLE" as the key.
    Return a GetSnakemakeRunDetails object with the filtered dataframe and the samplesheet filename.

    Parameters
    ----------
    inputs_obj : CLIparser
        The parsed command line arguments.

    Returns
    -------
    GetSnakemakeRunDetails
        The object containing the filtered samples dataframe and the samplesheet filename.
    """
    inputs_obj.samples_df = inputs_obj.samples_df[inputs_obj.samples_df["MATCH-REF"]]
    inputs_obj.samples_df.set_index("SAMPLE", inplace=True)

    inputs_obj.samples_dict = inputs_obj.samples_df.to_dict(orient="index")

    return GetSnakemakeRunDetails(inputs_obj, samplesheetfilename="samples_MatchRef")


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
    for i in zip(cols_left, cols_right):
        original_df[i[0]] = original_df.apply(
            lambda x: (
                override_df[i[1]][override_df["sample"] == x["SAMPLE"]].values[0]
                if x["SAMPLE"] in override_df["sample"].values and x[i[0]] != "NONE"
                else x[i[0]]
            ),
            axis=1,
        )
    return original_df


def process_match_ref(parsed_inputs: CLIparser) -> CLIparser:
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
    snakemakedetails = get_snakemake_run_details(inputs_obj_match_ref)
    log.info(
        f"{'='*20} [bold orange_red1] Starting Match-reference process [/bold orange_red1] {'='*20}"
    )

    if download_containers(snakemakedetails.snakemake_run_conf) != 0:
        log.error("Failed to download containers required for workflow.\nPlease check the logs and your settings for more information and try again later.")
        sys.exit(1)

    status = run_snakemake(inputs_obj_match_ref, snakemakedetails)

    workflow_state: Literal["Failed", "Success"] = (
        "Failed" if status is False else "Success"
    )
    if status is False:
        WriteReport(
            inputs_obj_match_ref.workdir,
            inputs_obj_match_ref.input_path,
            inputs_obj_match_ref.exec_start_path,
            inputs_obj_match_ref.user_config,
            snakemakedetails.snakemake_run_parameters,
            snakemakedetails.snakemake_run_conf,
            workflow_state,
        )
        exit(1)

    if snakemakedetails.snakemake_run_conf["dryrun"] is True:
        return parsed_inputs

    # status is True, this means that a file should exist at {workdir}/data/match_ref_results.pkl
    df = pd.read_pickle(f"{inputs_obj_match_ref.workdir}/data/match_ref_results.pkl")

    # ensure the "Primer_file" and "Feat_file" columns exist in df. If they don't, create them and set their values to "NONE"
    if "Primer_file" not in df.columns:
        df["Primer_file"] = "NONE"
    if "Feat_file" not in df.columns:
        df["Feat_file"] = "NONE"

    # keep only columns "sample", "Reference", "Reference_file", "Primer_file" and "Feat_file" from df
    filt_df = df[["sample", "Reference", "Reference_file", "Primer_file", "Feat_file"]]
    imploded_df = filt_df.groupby("sample").agg(set).reset_index()

    # replace all the sets with just one value to only that value (without the set) except for the Reference column
    imploded_df = replace_sets_to_singular_values(
        imploded_df, ["Reference_file", "Primer_file", "Feat_file"]
    )

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
