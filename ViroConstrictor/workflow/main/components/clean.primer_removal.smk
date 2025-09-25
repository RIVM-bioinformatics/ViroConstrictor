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
        f"{logdir}" "AmpliGone_{Virus}.{RefID}.{sample}.log",
    benchmark:
        f"{logdir}{bench}" "AmpliGone_{Virus}.{RefID}.{sample}.txt"
    threads: config["threads"]["PrimerRemoval"]
    resources:
        mem_mb=high_memory_job,
        runtime=55,
    params:
        amplicontype=config["amplicon_type"],
        primer_mismatch_rate=lambda wc: SAMPLES[wc.sample]["PRIMER-MISMATCH-RATE"],
        fragment_lookaround_size=lambda wc: f"--fragment-lookaround-size {int(size)}" if (size := SAMPLES[wc.sample].get("FRAGMENT-LOOKAROUND-SIZE")) else "",
        alignmentpreset=lambda wc: get_preset_parameter(
            preset_name=SAMPLES[wc.sample]["PRESET"],
            parameter_name=f"AmpliGone_AlignmentPreset_{config['platform']}",
        ),
        alignmentmatrix=lambda wc: get_preset_parameter(
            preset_name=SAMPLES[wc.sample]["PRESET"],
            parameter_name=f"AmpliGone_AlignmentMatrix_{config['platform']}",
        ),
        extrasettings=lambda wc: get_preset_parameter(
            preset_name=SAMPLES[wc.sample]["PRESET"],
            parameter_name=f"AmpliGone_ExtraSettings_{config['platform']}",
        ),
    shell:
        """
        echo "Running AmpliGone with primer file: {input.pr}" > {log}
        AmpliGone \
            -i {input.fq} \
            -ref {input.ref} -pr {input.pr} \
            -o {output.fq} \
            -at {params.amplicontype} \
            --error-rate {params.primer_mismatch_rate} \
            {params.fragment_lookaround_size} \
            --export-primers {output.ep} \
            {params.alignmentpreset} {params.alignmentmatrix} {params.extrasettings} \
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
        runtime=55,
    shell:
        """
        cp {input} {output.fq}
        """
