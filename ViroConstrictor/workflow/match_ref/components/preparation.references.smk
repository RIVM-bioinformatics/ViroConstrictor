
rule filter_references:
    input:
        lambda wc: SAMPLES[wc.sample]["REFERENCE"],
    output:
        temp(f"{datadir}{matchref}{wc_folder}" "{sample}_refs.fasta"),
    resources:
        mem_mb=low_memory_job,
        runtime=low_runtime_job,
    threads: 1
    log:
        f"{logdir}prepare_refs" "{Virus}.{segment}.{sample}.log",
    benchmark:
        f"{logdir}{bench}MR_prepare_refs" "{Virus}.{segment}.{sample}.txt"
    conda:
        workflow_environment_path("mr_scripts.yaml")
    container:
        f"{container_base_path}/viroconstrictor_mr_scripts_{get_hash('mr_scripts')}.sif"
    params:
        script="-m match_ref.scripts.filter_references",
        pythonpath=f'{Path(workflow.basedir).parent}'
    shell:
        """
        PYTHONPATH={params.pythonpath} \
        python {params.script} \
        --input {input} \
        --output {output} \
        --wildcard_segment {wildcards.segment} >> {log} 2>&1
        """
