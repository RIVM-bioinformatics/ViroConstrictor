
rule vcf_to_tsv:
    input:
        vcf=rules.trueconsense.output.vcf,
    output:
        tsv=temp(f"{datadir}{wc_folder}{aln}{vf}" "{sample}.tsv"),
    conda:
        workflow_environment_path("core_scripts.yaml")
    container:
        f"{container_base_path}/viroconstrictor_core_scripts_{get_hash('core_scripts')}.sif"
    threads: config["threads"]["Index"]
    resources:
        mem_mb=low_memory_job,
        runtime=55,
    log:
        f"{logdir}" "vcf_to_tsv_{Virus}.{RefID}.{sample}.log",
    params:
        script="-m main.scripts.vcf_to_tsv",
        pythonpath = f'{Path(workflow.basedir).parent}',
    shell:
        """
        PYTHONPATH={params.pythonpath} \
        python {params.script} \
        --input {input.vcf} \
        --output {output.tsv} \
        --sample {wildcards.sample} >> {log} 2>&1
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
        workflow_environment_path("core_scripts.yaml")
    container:
        f"{container_base_path}/viroconstrictor_core_scripts_{get_hash('core_scripts')}.sif"
    params:
        script="-m main.scripts.boc",
        pythonpath = f'{Path(workflow.basedir).parent}',
    shell:
        """
        PYTHONPATH={params.pythonpath} \
        python {params.script} \
        --input {input.reference} \
        --sample {wildcards.sample} \
        --coverage {input.coverage} \
        --output {output}
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
        workflow_environment_path("core_scripts.yaml")
    container:
        f"{container_base_path}/viroconstrictor_core_scripts_{get_hash('core_scripts')}.sif"
    params:
        script="-m main.scripts.amplicon_covs",
        pythonpath = f'{Path(workflow.basedir).parent}',
    shell:
        """
        PYTHONPATH={params.pythonpath} \
        python {params.script} \
        --input {input.pr} \
        --coverages {input.cov} \
        --key {wildcards.sample} \
        --output {output} > {log} 2>&1
        """
