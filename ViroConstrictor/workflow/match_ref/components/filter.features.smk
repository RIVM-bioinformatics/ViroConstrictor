
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
        runtime=low_runtime_job,
    log:
        f"{logdir}FilterGFF_" "{sample}.log",
    benchmark:
        f"{logdir}{bench}FilterGFF_" "{sample}.txt",
    conda:
        workflow_environment_path("mr_scripts.yaml")
    container:
        f"{container_base_path}/viroconstrictor_mr_scripts_{get_hash('mr_scripts')}.sif"
    params:
        script="-m match_ref.scripts.filter_gff",
        pythonpath=f'{Path(workflow.basedir).parent}'
    shell:
        """
        PYTHONPATH={params.pythonpath} \
        python {params.script} \
        --refdata {input.refdata} \
        --input {input.gff} \
        --output {output.gff} \
        --updatedstats {output.groupedstats} >> {log} 2>&1
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
        runtime=low_runtime_job,
    shell:
        """
        cp {input.refdata} {output.groupedstats}
        touch {output.gff}
        """
