
rule filter_best_matching_ref:
    input:
        stats=rules.count_mapped_reads.output,
        ref=rules.filter_references.output,
    output:
        filtref=temp(f"{datadir}{matchref}{wc_folder}" "{sample}_best_ref.fasta"),
        filtcount=temp(f"{datadir}{matchref}{wc_folder}" "{sample}_best_ref.csv"),
    conda:
        workflow_environment_path("Scripts.yaml")
    container:
        f"{container_base_path}/viroconstrictor_scripts_{get_hash('Scripts')}.sif"
    threads: 1
    resources:
        mem_mb=low_memory_job,
        runtime=55,
    log:
        f"{logdir}FilterBR_" "{Virus}.{segment}.{sample}.log",
    params:
        script=(
            workflow_script_path("scripts/filter_best_matching_ref.py")
            if (
                DeploymentMethod.CONDA
                in workflow.deployment_settings.deployment_method
            )
            is True
            else "/match_ref_scripts/filter_best_matching_ref.py"
        ),
    shell:
        """
        python {params.script} {input.stats} {input.ref} {output.filtref} {output.filtcount} >> {log} 2>&1
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
        workflow_environment_path("Scripts.yaml")
    container:
        f"{container_base_path}/viroconstrictor_scripts_{get_hash('Scripts')}.sif"
    threads: 1
    resources:
        mem_mb=low_memory_job,
        runtime=55,
    log:
        f"{logdir}GroupRefs_" "{sample}.log",
    params:
        script=(
            workflow_script_path("scripts/group_refs.py")
            if (
                DeploymentMethod.CONDA
                in workflow.deployment_settings.deployment_method
            )
            is True
            else "/match_ref_scripts/group_refs.py"
        ),
    shell:
        """
        python {params.script} "{input.ref}" "{input.stats}" {output.groupedrefs} {output.groupedstats} {wildcards.sample} >> {log} 2>&1
        """
