def get_primers_output(wildcards):
    sample_primers = SAMPLES[wildcards.sample]["PRIMERS"]
    if sample_primers == "NONE":
        return ""
    elif sample_primers.endswith(".bed"):
        return rules.filter_primer_bed.output.bed
    else:
        return rules.prepare_primers.output.bed


rule ampligone:
    input:
        fq=rules.qc_filter.output.fq,
        pr=lambda wc: get_primers_output(wc),
        ref=rules.prepare_refs.output,
    output:
        fq=f"{datadir}{wc_folder}{cln}{prdir}" "{sample}.fastq",
        ep=f"{datadir}{wc_folder}{prim}" "{sample}_removedprimers.bed",
    conda:
        workflow_environment_path("Clean.yaml")
    container:
        f"{container_base_path}/viroconstrictor_clean_{get_hash('Clean')}.sif"
    log:
        f"{logdir}{get_rule_name()}/" "AmpliGone_{Virus}.{RefID}.{sample}.log",
    threads: config["threads"]["PrimerRemoval"]
    resources:
        mem_mb=high_memory_job,
        runtime=high_runtime_job,
    params:
        binary = lambda wc: get_binary(
            stage_identifier=VC_STAGE,
            preset_name=SAMPLES[wc.sample]["PRESET"],
            task="primer_removal",
            platform=config["platform"],
        ),
        flags = lambda wc: get_flags(
            stage_identifier=VC_STAGE,
            preset_name=SAMPLES[wc.sample]["PRESET"],
            task="primer_removal",
            platform=config["platform"],
        ),
    shell:
        """
        echo "Running AmpliGone with primer file: {input.pr}" > {log}
        {params.binary} \
            -i {input.fq} \
            -ref {input.ref} -pr {input.pr} \
            -o {output.fq} \
            --export-primers {output.ep} \
            {params.flags} \
            -to \
            -t {threads} >> {log} 2>&1
        """


ruleorder: ampligone > move_fastq


# If no primers are given (e.g. with illumina runs), this rule makes sure the fastq's end up in the right place
rule move_fastq:
    input:
        rules.qc_filter.output.fq,
    output:
        fq=f"{datadir}{wc_folder}{cln}{prdir}" "{sample}.fastq",
        ep=touch(f"{datadir}{wc_folder}{prim}" "{sample}_removedprimers.bed"),
    resources:
        mem_mb=low_memory_job,
        runtime=low_runtime_job,
    shell:
        """
        cp {input} {output.fq}
        """
