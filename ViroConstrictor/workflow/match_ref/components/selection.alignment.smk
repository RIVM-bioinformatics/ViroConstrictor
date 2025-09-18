if config["platform"] in ["nanopore", "iontorrent"] or (
    config["platform"] == "illumina" and config["unidirectional"] is True
):
    base_mm2_preset = (
        "-ax sr" if config["platform"] in ["iontorrent", "illumina"] else "-ax map-ont"
    )

    rule align_to_refs:
        input:
            ref=rules.filter_references.output,
            fq=lambda wc: SAMPLES[wc.sample]["INPUTFILE"],
        output:
            bam=temp(f"{datadir}{matchref}{wc_folder}" "{sample}.bam"),
            index=temp(f"{datadir}{matchref}{wc_folder}" "{sample}.bam.bai"),
        conda:
            workflow_environment_path("Alignment.yaml")
        container:
            f"{container_base_path}/viroconstrictor_alignment_{get_hash('Alignment')}.sif"
        log:
            f"{logdir}AlignMR_" "{Virus}.{segment}.{sample}.log",
        benchmark:
            f"{logdir}{bench}AlignMR_" "{Virus}.{segment}.{sample}.txt"
        threads: config["threads"]["Alignments"]
        resources:
            mem_mb=medium_memory_job,
            runtime=55,
        params:
            mapthreads=config["threads"]["Alignments"] - 1,
            mm2_alignment_preset=base_mm2_preset,
            minimap2_base_setting=lambda wc: get_preset_parameter(
                preset_name=SAMPLES[wc.sample]["PRESET"],
                parameter_name="Minimap2_Settings_Base",
            ),
            minimap2_extra_setting=lambda wc: get_preset_parameter(
                preset_name=SAMPLES[wc.sample]["PRESET"],
                parameter_name="Minimap2_Settings",
            ),
            minimap2_alignmentparams=lambda wc: get_preset_parameter(
                preset_name=SAMPLES[wc.sample]["PRESET"],
                parameter_name=f"Minimap2_AlignmentParams_{config['platform']}",
            ),
            samtools_standard_filters=lambda wc: get_preset_parameter(
                preset_name=SAMPLES[wc.sample]["PRESET"],
                parameter_name="Samtools_Filters_Base",
            ),
            samtools_extra_filters=lambda wc: get_preset_parameter(
                preset_name=SAMPLES[wc.sample]["PRESET"],
                parameter_name=f"Samtools_Filters_{config['platform']}",
            ),
        shell:
            """
            minimap2 {params.mm2_alignment_preset} {params.minimap2_base_setting} {params.minimap2_extra_setting} {params.minimap2_alignmentparams} -t {params.mapthreads} {input.ref} {input.fq} 2>> {log} |\
            samtools view -@ {threads} {params.samtools_standard_filters} {params.samtools_extra_filters} -uS 2>> {log} |\
            samtools sort -o {output.bam} >> {log} 2>&1
            samtools index {output.bam} >> {log} 2>&1
            """


if config["platform"] == "illumina" and config["unidirectional"] is False:
    base_mm2_preset = "-ax sr"

    rule align_to_refs:
        input:
            ref=rules.filter_references.output,
            fq1=lambda wildcards: SAMPLES[wildcards.sample]["R1"],
            fq2=lambda wildcards: SAMPLES[wildcards.sample]["R2"],
        output:
            bam=temp(f"{datadir}{matchref}{wc_folder}" "{sample}.bam"),
            index=temp(f"{datadir}{matchref}{wc_folder}" "{sample}.bam.bai"),
        conda:
            workflow_environment_path("Alignment.yaml")
        container:
            f"{container_base_path}/viroconstrictor_alignment_{get_hash('Alignment')}.sif"
        log:
            f"{logdir}" "AlignMR_{Virus}.{segment}.{sample}.log",
        benchmark:
            f"{logdir}{bench}" "AlignMR_{Virus}.{segment}.{sample}.log"
        threads: config["threads"]["Alignments"]
        resources:
            mem_mb=medium_memory_job,
            runtime=55,
        params:
            mapthreads=config["threads"]["Alignments"] - 1,
            mm2_alignment_preset=base_mm2_preset,
            minimap2_base_setting=lambda wc: get_preset_parameter(
                preset_name=SAMPLES[wc.sample]["PRESET"],
                parameter_name="Minimap2_Settings_Base",
            ),
            minimap2_extra_setting=lambda wc: get_preset_parameter(
                preset_name=SAMPLES[wc.sample]["PRESET"],
                parameter_name="Minimap2_Settings",
            ),
            minimap2_alignmentparams=lambda wc: get_preset_parameter(
                preset_name=SAMPLES[wc.sample]["PRESET"],
                parameter_name=f"Minimap2_AlignmentParams_{config['platform']}",
            ),
            samtools_standard_filters=lambda wc: get_preset_parameter(
                preset_name=SAMPLES[wc.sample]["PRESET"],
                parameter_name="Samtools_Filters_Base",
            ),
            samtools_extra_filters=lambda wc: get_preset_parameter(
                preset_name=SAMPLES[wc.sample]["PRESET"],
                parameter_name=f"Samtools_Filters_{config['platform']}",
            ),
        shell:
            """
            minimap2 {params.mm2_alignment_preset} {params.minimap2_base_setting} {params.minimap2_extra_setting} {params.minimap2_alignmentparams} -t {params.mapthreads} {input.ref} {input.fq1} {input.fq2} 2>> {log} |\
            samtools view -@ {threads} {params.samtools_standard_filters} {params.samtools_extra_filters} -uS 2>> {log} |\
            samtools sort -o {output.bam} >> {log} 2>&1
            samtools index {output.bam} >> {log} 2>&1
            """


rule count_mapped_reads:
    input:
        bam=rules.align_to_refs.output.bam,
        index=rules.align_to_refs.output.index,
    output:
        temp(f"{datadir}{matchref}{wc_folder}" "{sample}_count.csv"),
    conda:
        workflow_environment_path("mr_scripts.yaml")
    container:
        f"{container_base_path}/viroconstrictor_mr_scripts_{get_hash('mr_scripts')}.sif"
    threads: 1
    resources:
        mem_mb=low_memory_job,
        runtime=55,
    log:
        f"{logdir}CountMR_" "{Virus}.{segment}.{sample}.log",
    params:
        script="-m match_ref.scripts.count_mapped_reads",
        pythonpath=f'{Path(workflow.basedir).parent}'
    shell:
        """
        PYTHONPATH={params.pythonpath} \
        python {params.script} \
        --input {input.bam} \
        --output {output} >> {log} 2>&1
        """
