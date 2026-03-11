
rule concat_sequences:
    input:  # we don't use group_items here as mincov is in the filename
        lambda wc: (
            f"{datadir}{wc_folder}{cons}{seqs}{sample}.fa"
            for sample in p_space.dataframe.loc[
                (p_space.Virus == wc.Virus) & (p_space.RefID == wc.RefID)
            ]["sample"]
        ),
    output:
        f"{res}{wc_folder}consensus.fasta",
    resources:
        mem_mb=low_memory_job,
        runtime=low_runtime_job,
    shell:
        "cat {input} >> {output}"


def group_aminoacids_inputs(wildcards):
    filtered_df = samples_df.loc[samples_df["AA_FEAT_NAMES"].notnull()]
    filtered_vir_list = list(filtered_df["Virus"].unique())

    struct = {}
    for i in filtered_vir_list:
        select_samples = list(
            samples_df.loc[samples_df["Virus"] == i]["sample"].unique()
        )
        # for x in select_samples:
        #     y = samples_df.loc[(samples_df["Virus"] == i) & (samples_df["sample"] == x)]["RefID"].unique()
        #     print(y)
        # select_refIDs = list(samples_df.loc[samples_df["Virus"] == i]["RefID"].unique())
        # print(select_refIDs)

        # create a dictionary of dictionaries for each virus, with 'i' as the primary key and sample as the secondary key having a list of refIDs as the value
        struct[i] = {
            sample: list(
                samples_df.loc[
                    (samples_df["Virus"] == i) & (samples_df["sample"] == sample)
                ]["RefID"].unique()
            )
            for sample in select_samples
        }
    file_list = []
    for virus, sample in struct.items():
        for sample, refid in sample.items():
            for ref in refid:
                file_list.append(
                    f"{datadir}Virus~{virus}/RefID~{ref}/{amino}{sample}/aa.faa"
                )
    return file_list


# this rule cannot and should not run in a separate environment/container as the sole purpose is to transfer data of the paramspace into something that can then be used in a later rule.
# TODO: this rule should probably be moved to core workflow.smk as it is not specific to concatenations.
rule make_pickle:
    output:
        temp(f"{datadir}sampleinfo.pkl"),
    resources:
        mem_mb=low_memory_job,
        runtime=low_runtime_job,
    threads: 1
    run:
        import pandas as pd
        # select rows with aminoacid feature names and write pickle in the current (Snakemake) runtime
        df = samples_df[~samples_df["AA_FEAT_NAMES"].isnull()]
        df.to_pickle(output[0], compression=None)


rule concat_aminoacids:
    input:
        files=lambda wildcards: group_aminoacids_inputs(wildcards),
        sampleinfo=rules.make_pickle.output,
    output:
        list_aminoacid_result_outputs(samples_df),
    log:
        f"{logdir}concat_aminoacids.log",
    resources:
        mem_mb=low_memory_job,
        runtime=low_runtime_job,
    conda:
        workflow_environment_path("core_scripts.yaml")
    container:
        f"{container_base_path}/viroconstrictor_core_scripts_{get_hash('core_scripts')}.sif"
    threads: 1
    params:
        script="-m main.scripts.group_aminoacids",
        pythonpath=f'{Path(workflow.basedir).parent}',
    shell:
        """
        PYTHONPATH={params.pythonpath} \
        python {params.script} \
        --input "{input.files}" \
        --output "{output}" \
        --space {input.sampleinfo} > {log} 2>&1
        """

def group_items(wildcards, folder, filename):
    filtered_virus = p_space.dataframe.loc[
        p_space.dataframe["Virus"] == wildcards.Virus
    ]
    filtered_refid = filtered_virus.loc[filtered_virus["RefID"] == wildcards.RefID]
    return [f"{folder}{item}{filename}" for item in list(filtered_refid["sample"])]


rule concat_tsv_coverages:
    input:
        lambda wildcards: group_items(
            wildcards, folder=f"{datadir}{wc_folder}{aln}{vf}", filename=".tsv"
        ),
    output:
        f"{res}{wc_folder}mutations.tsv",
    resources:
        mem_mb=low_memory_job,
        runtime=low_runtime_job,
    shell:
        """
        echo -e 'Sample\tReference_ID\tPosition\tReference_Base\tVariant_Base\tDepth' > {output}
        cat {input} >> {output}
        """


rule concat_boc:
    input:
        lambda wildcards: group_items(
            wildcards, folder=f"{datadir}{wc_folder}{boc}", filename=".tsv"
        ),
    output:
        f"{res}{wc_folder}Width_of_coverage.tsv",
    resources:
        mem_mb=low_memory_job,
        runtime=low_runtime_job,
    shell:
        """
        echo -e "Sample_name\tWidth_at_mincov_1\tWidth_at_mincov_5\tWidth_at_mincov_10\tWidth_at_mincov_50\tWidth_at_mincov_100" > {output}
        cat {input} >> {output}
        """


rule concat_amplicon_cov:
    input:
        pr = lambda wildcards: group_items(
            wildcards,
            folder=f"{datadir}{wc_folder}{prim}",
            filename="_ampliconcoverage.csv",
        ),
    output:
        f"{res}{wc_folder}Amplicon_coverage.csv",
    resources:
        mem_mb=low_memory_job,
        runtime=low_runtime_job,
    conda:
        workflow_environment_path("core_scripts.yaml")
    container:
        f"{container_base_path}/viroconstrictor_core_scripts_{get_hash('core_scripts')}.sif"
    log:
        f"{logdir}concat_amplicon_cov_" "{Virus}.{RefID}.log",
    threads: 1
    params:
        script="-m main.scripts.concat_amplicon_covs",
        pythonpath=f'{Path(workflow.basedir).parent}',
    shell:
        """
        PYTHONPATH={params.pythonpath} \
        python {params.script} \
        --input "empty" \
        --input_coverages {input.pr} \
        --output {output}
        """
