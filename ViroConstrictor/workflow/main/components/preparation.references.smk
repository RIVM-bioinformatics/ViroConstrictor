

rule prepare_refs:
    input:
        lambda wc: SAMPLES[wc.sample]["REFERENCE"],
    output:
        f"{datadir}{wc_folder}" "{sample}_reference.fasta",
    resources:
        mem_mb=low_memory_job,
        runtime=55,
    threads: 1
    log:
        f"{logdir}prepare_refs_" "{Virus}.{RefID}.{sample}.log",
    benchmark:
        f"{logdir}{bench}prepare_refs_" "{Virus}.{RefID}.{sample}.txt"
    conda:
        workflow_environment_path("Scripts.yaml")
    container:
        f"{container_base_path}/viroconstrictor_scripts_{get_hash('Scripts')}.sif"
    params:
        script=(
            workflow_script_path("scripts/prepare_refs.py")
            if (
                DeploymentMethod.CONDA
                in workflow.deployment_settings.deployment_method
            )
            is True
            else "/scripts/prepare_refs.py"
        ),
    shell:
        """
        python {params.script} {input} {output} {wildcards.RefID} > {log}
        """
