

rule prepare_refs:
    input:
        lambda wc: SAMPLES[wc.sample]["REFERENCE"],
    output:
        f"{datadir}{wc_folder}" "{sample}_reference.fasta",
    resources:
        mem_mb=low_memory_job,
        runtime=low_runtime_job,
    threads: 1
    log:
        f"{logdir}prepare_refs_" "{Virus}.{RefID}.{sample}.log",
    benchmark:
        f"{logdir}{bench}prepare_refs_" "{Virus}.{RefID}.{sample}.txt"
    conda:
        workflow_environment_path("core_scripts.yaml")
    container:
        f"{container_base_path}/viroconstrictor_core_scripts_{get_hash('core_scripts')}.sif"
    params:
        script="-m main.scripts.prepare_refs",
        pythonpath=f'{Path(workflow.basedir).parent}'
    shell:
        """
        PYTHONPATH={params.pythonpath} \
        python {params.script} \
        --input {input} \
        --output {output} \
        --reference_id {wildcards.RefID} > {log}
        """