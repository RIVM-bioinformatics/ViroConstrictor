
rule vcf_to_tsv:
    input:
        vcf=rules.trueconsense.output.vcf,
    output:
        tsv=temp(f"{datadir}{wc_folder}{aln}{vf}" "{sample}.tsv"),
    conda:
        workflow_environment_path("Scripts.yaml")
    container:
        f"{container_base_path}/viroconstrictor_scripts_{get_hash('Scripts')}.sif"
    threads: config["threads"]["Index"]
    resources:
        mem_mb=low_memory_job,
        runtime=55,
    log:
        f"{logdir}" "vcf_to_tsv_{Virus}.{RefID}.{sample}.log",
    params:
        script=(
            workflow_script_path("scripts/vcf_to_tsv.py")
            if (
                DeploymentMethod.CONDA
                in workflow.deployment_settings.deployment_method
            )
            is True
            else "/scripts/vcf_to_tsv.py"
        ),
    shell:
        """
        python {params.script} {input.vcf} {output.tsv} {wildcards.sample} >> {log} 2>&1
        """


rule get_breadth_of_coverage:
    input:
        reference=rules.prepare_refs.output,
        coverage=rules.trueconsense.output.cov,
    output:
        temp(f"{datadir}{wc_folder}{boc}" "{sample}.tsv"),
    resources:
        mem_mb=low_memory_job,
        runtime=55,
    conda:
        workflow_environment_path("Scripts.yaml")
    container:
        f"{container_base_path}/viroconstrictor_scripts_{get_hash('Scripts')}.sif"
    params:
        script=(
            workflow_script_path("scripts/boc.py")
            if (
                DeploymentMethod.CONDA
                in workflow.deployment_settings.deployment_method
            )
            is True
            else "/scripts/boc.py"
        ),
    shell:
        """
        python {params.script} {input.reference} {wildcards.sample} {input.coverage} {output}
        """


rule calculate_amplicon_cov:
    input:
        pr=f"{datadir}{wc_folder}{prim}" "{sample}_primers.bed",
        cov=rules.trueconsense.output.cov,
    output:
        f"{datadir}{wc_folder}{prim}" "{sample}_ampliconcoverage.csv",
    log:
        f"{logdir}" "calculate_amplicon_cov_{Virus}.{RefID}.{sample}.log",
    benchmark:
        f"{logdir}{bench}" "calculate_amplicon_cov_{Virus}.{RefID}.{sample}.txt"
    resources:
        mem_mb=low_memory_job,
        runtime=55,
    conda:
        workflow_environment_path("Scripts.yaml")
    container:
        f"{container_base_path}/viroconstrictor_scripts_{get_hash('Scripts')}.sif"
    params:
        script=(
            workflow_script_path("scripts/amplicon_covs.py")
            if (
                DeploymentMethod.CONDA
                in workflow.deployment_settings.deployment_method
            )
            is True
            else "/scripts/amplicon_covs.py"
        ),
    shell:
        """
        python {params.script} \
        --primers {input.pr} \
        --coverages {input.cov} \
        --key {wildcards.sample} \
        --output {output} > {log} 2>&1
        """
