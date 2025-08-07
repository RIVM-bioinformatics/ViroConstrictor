
rule filter_references:
    input:
        lambda wc: SAMPLES[wc.sample]["REFERENCE"],
    output:
        temp(f"{datadir}{matchref}{wc_folder}" "{sample}_refs.fasta"),
    resources:
        mem_mb=low_memory_job,
        runtime=55,  # TODO: this just sets a default runtime of 55 minutes. However this should be dynamic just like the memory.
    threads: 1
    log:
        f"{logdir}prepare_refs" "{Virus}.{segment}.{sample}.log",
    benchmark:
        f"{logdir}{bench}MR_prepare_refs" "{Virus}.{segment}.{sample}.txt"
    conda:
        workflow_environment_path("Scripts.yaml")
    container:
        f"{container_base_path}/viroconstrictor_scripts_{get_hash('Scripts')}.sif"
    params:
        script=(
            workflow_script_path("scripts/filter_references.py")
            if (
                DeploymentMethod.CONDA
                in workflow.deployment_settings.deployment_method
            )
            is True
            else "/match_ref_scripts/filter_references.py"
        ),
    shell:
        """
        python {params.script} {input} {output} {wildcards.segment} >> {log} 2>&1
        """
