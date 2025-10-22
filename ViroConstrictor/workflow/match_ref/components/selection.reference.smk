
rule filter_best_matching_ref:
    input:
        stats=rules.count_mapped_reads.output,
        ref=rules.filter_references.output,
    output:
        filtref=temp(f"{datadir}{matchref}{wc_folder}" "{sample}_best_ref.fasta"),
        filtcount=temp(f"{datadir}{matchref}{wc_folder}" "{sample}_best_ref.csv"),
    conda:
        workflow_environment_path("mr_scripts.yaml")
    container:
        f"{container_base_path}/viroconstrictor_mr_scripts_{get_hash('mr_scripts')}.sif"
    threads: 1
    resources:
        mem_mb=low_memory_job,
        runtime=55,
    log:
        f"{logdir}FilterBR_" "{Virus}.{segment}.{sample}.log",
    params:
        script="-m match_ref.scripts.filter_best_matching_ref",
        pythonpath=f'{Path(workflow.basedir).parent}'
    shell:
        """
        PYTHONPATH={params.pythonpath} \
        python {params.script} \
        --input {input.stats} \
        --inputref {input.ref} \
        --filtref {output.filtref} \
        --output {output.filtcount} >> {log} 2>&1
        """


rule group_and_rename_refs:
    input:
        ref=lambda wildcards: expand(
            f"{datadir}{matchref}{wc_folder}" "{sample}_best_ref.fasta",
            zip,
            Virus=p_space.dataframe.loc[
                p_space.dataframe["sample"] == wildcards.sample, "Virus"
            ],
            segment=p_space.dataframe.loc[
                p_space.dataframe["sample"] == wildcards.sample, "segment"
            ],
            allow_missing=True,
        ),
        stats=lambda wildcards: expand(
            f"{datadir}{matchref}{wc_folder}" "{sample}_best_ref.csv",
            zip,
            Virus=p_space.dataframe.loc[
                p_space.dataframe["sample"] == wildcards.sample, "Virus"
            ],
            segment=p_space.dataframe.loc[
                p_space.dataframe["sample"] == wildcards.sample, "segment"
            ],
            allow_missing=True,
        ),
    output:
        groupedrefs=f"{datadir}{matchref}" "{sample}_refs.fasta",
        groupedstats=temp(f"{datadir}{matchref}" "{sample}_refs.csv"),
    conda:
        workflow_environment_path("mr_scripts.yaml")
    container:
        f"{container_base_path}/viroconstrictor_mr_scripts_{get_hash('mr_scripts')}.sif"
    threads: 1
    resources:
        mem_mb=low_memory_job,
        runtime=55,
    log:
        f"{logdir}GroupRefs_" "{sample}.log",
    params:
        script="-m match_ref.scripts.group_refs",
        pythonpath=f'{Path(workflow.basedir).parent}'
    shell:
        """
        PYTHONPATH={params.pythonpath} \
        python {params.script} \
        --input "empty" \
        --input_refs {input.ref} \
        --input_stats {input.stats} \
        --output {output.groupedrefs} \
        --output_stats {output.groupedstats} \
        --sample {wildcards.sample} >> {log} 2>&1
        """
