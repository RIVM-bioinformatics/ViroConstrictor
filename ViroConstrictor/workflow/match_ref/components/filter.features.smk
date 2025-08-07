
rule filter_gff:
    input:
        refdata=rules.group_and_rename_refs.output.groupedstats,
        gff=lambda wc: file if (file := SAMPLES[wc.sample]["FEATURES"]) else "",
    output:
        gff=f"{datadir}{matchref}" "{sample}_feats.gff",
        groupedstats=temp(f"{datadir}{matchref}" "{sample}_step1.csv"),
    threads: 1
    resources:
        mem_mb=low_memory_job,
        runtime=55,
    log:
        f"{logdir}FilterGFF_" "{sample}.log",
    conda:
        workflow_environment_path("Scripts.yaml")
    container:
        f"{container_base_path}/viroconstrictor_scripts_{get_hash('Scripts')}.sif"
    params:
        script=(
            workflow_script_path("scripts/filter_gff.py")
            if (
                DeploymentMethod.CONDA
                in workflow.deployment_settings.deployment_method
            )
            is True
            else "/match_ref_scripts/filter_gff.py"
        ),
    shell:
        """
        python {params.script} {input.refdata} {input.gff} {output.gff} {output.groupedstats} >> {log} 2>&1
        """


ruleorder: filter_gff > touch_gff


rule touch_gff:
    input:
        refdata=rules.group_and_rename_refs.output.groupedstats,
    output:
        gff=f"{datadir}{matchref}" "{sample}_feats.gff",
        groupedstats=temp(f"{datadir}{matchref}" "{sample}_step1.csv"),
    threads: 1
    resources:
        mem_mb=low_memory_job,
        runtime=55,
    shell:
        """
        cp {input.refdata} {output.groupedstats}
        touch {output.gff}
        """
