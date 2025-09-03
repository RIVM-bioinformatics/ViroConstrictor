
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
        runtime=55,
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
rule make_pickle:
    output:
        temp(f"{datadir}sampleinfo.pkl"),
    resources:
        mem_mb=low_memory_job,
        runtime=55,
    threads: 1
    params:
        space=lambda wc: __import__("codecs")
        .encode(
            __import__("pickle").dumps(
                samples_df[~samples_df["AA_FEAT_NAMES"].isnull()]
            ),
            "base64",
        )
        .decode(),
    shell:
        """
        python -c "import pickle, codecs, pandas; df = pickle.loads(codecs.decode('''{params.space}'''.encode(), 'base64')); df.to_pickle('{output}', compression=None)"
        """


rule concat_aminoacids:
    input:
        files=lambda wildcards: group_aminoacids_inputs(wildcards),
        sampleinfo=rules.make_pickle.output,
    output:
        list_aminoacid_result_outputs(samples_df),
    resources:
        mem_mb=low_memory_job,
        runtime=55,
    conda:
        workflow_environment_path("Scripts.yaml")
    container:
        f"{container_base_path}/viroconstrictor_scripts_{get_hash('Scripts')}.sif"
    threads: 1
    params:
        script="-m scripts.group_aminoacids",
    shell:
        """
        PYTHONPATH={workflow.basedir} \
        python {params.script} \
        --input "{input.files}" \
        --output "{output}" \
        --space {input.sampleinfo}
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
        runtime=55,
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
        runtime=55,
    shell:
        """
        echo -e "Sample_name\tWidth_at_mincov_1\tWidth_at_mincov_5\tWidth_at_mincov_10\tWidth_at_mincov_50\tWidth_at_mincov_100" > {output}
        cat {input} >> {output}
        """


rule concat_amplicon_cov:
    input:
        lambda wildcards: group_items(
            wildcards,
            folder=f"{datadir}{wc_folder}{prim}",
            filename="_ampliconcoverage.csv",
        ),
    output:
        f"{res}{wc_folder}Amplicon_coverage.csv",
    resources:
        mem_mb=low_memory_job,
        runtime=55,
    conda:
        workflow_environment_path("Scripts.yaml")
    container:
        f"{container_base_path}/viroconstrictor_scripts_{get_hash('Scripts')}.sif"
    params:
        script = "-m scripts.amplicon_covs",
    shell:
        """
        PYTHONPATH={workflow.basedir} \
        python {params.script} \
        --input {input.pr} \
        --coverages {input.cov} \
        --key {wildcards.sample} \
        --output {output} > {log} 2>&1
        """
