
rule prepare_primers:
    input:
        prm=lambda wc: (
            ""
            if SAMPLES[wc.sample]["PRIMERS"].endswith(".bed")
            else SAMPLES[wc.sample]["PRIMERS"]
        ),
        ref=rules.prepare_refs.output,
    output:
        bed=f"{datadir}{wc_folder}{prim}" "{sample}_primers.bed",
    resources:
        mem_mb=low_memory_job,
        runtime=55,
    log:
        f"{logdir}prepare_primers_" "{Virus}.{RefID}.{sample}.log",
    benchmark:
        f"{logdir}{bench}prepare_primers_" "{Virus}.{RefID}.{sample}.txt"
    params:
        pr_mm_rate=lambda wc: SAMPLES[wc.sample]["PRIMER-MISMATCH-RATE"],
    conda:
        workflow_environment_path("Clean.yaml")
    container:
        f"{container_base_path}/viroconstrictor_clean_{get_hash('Clean')}.sif"
    shell:
        """
        python -m AmpliGone.fasta2bed \
            --primers {input.prm} \
            --reference {input.ref} \
            --output {output.bed} \
            --primer-mismatch-rate {params.pr_mm_rate} \
            --verbose > {log}
        """


ruleorder: prepare_primers > filter_primer_bed > create_empty_primers


rule filter_primer_bed:
    input:
        prm=lambda wc: SAMPLES[wc.sample]["PRIMERS"],
    output:
        bed=f"{datadir}{wc_folder}{prim}" "{sample}_primers.bed",
    resources:
        mem_mb=low_memory_job,
        runtime=55,
    log:
        f"{logdir}prepare_primers_" "{Virus}.{RefID}.{sample}.log",
    benchmark:
        f"{logdir}{bench}prepare_primers_" "{Virus}.{RefID}.{sample}.txt"
    params:
        script=(
            workflow_script_path("scripts/filter_bed_input.py")
            if (
                DeploymentMethod.CONDA
                in workflow.deployment_settings.deployment_method
            )
            is True
            else "/scripts/filter_bed_input.py"
        ),
    conda:
        workflow_environment_path("Scripts.yaml")
    container:
        f"{container_base_path}/viroconstrictor_scripts_{get_hash('Scripts')}.sif"
    shell:
        """
        python {params.script} {input.prm} {output.bed} {wildcards.RefID}
        """


rule create_empty_primers:
    output:
        bed=touch(f"{datadir}{wc_folder}{prim}" "{sample}_primers.bed"),
    resources:
        mem_mb=low_memory_job,
    log:
        f"{logdir}prepare_primers_" "{Virus}.{RefID}.{sample}.log",
    shell:
        """
        echo "Created empty primer bed file for NONE primers" > {log}
        """
