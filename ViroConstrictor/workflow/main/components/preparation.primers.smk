
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
    conda:
        workflow_environment_path("core_scripts.yaml")
    container:
        f"{container_base_path}/viroconstrictor_core_scripts_{get_hash('core_scripts')}.sif"
    params:
        pythonpath = f'{Path(workflow.basedir).parent}',
        script="-m main.scripts.filter_bed_input",
    shell:
        """
        PYTHONPATH={params.pythonpath} \
        python {params.script} \
        --input {input.prm} \
        --output {output.bed} \
        --reference_id {wildcards.RefID}
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
