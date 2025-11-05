if config["platform"] in ["nanopore", "iontorrent"] or (
    config["platform"] == "illumina" and config["unidirectional"] is True
):
    base_mm2_preset = (
        "-ax sr" if config["platform"] in ["iontorrent", "illumina"] else "-ax map-ont"
    )

    rule remove_adapters_p1:
        input:
            ref=rules.prepare_refs.output,
            fq=lambda wc: SAMPLES[wc.sample]["INPUTFILE"],
        output:
            bam=f"{datadir}{wc_folder}{cln}{raln}" "{sample}.bam",
            index=f"{datadir}{wc_folder}{cln}{raln}" "{sample}.bam.bai",
        conda:
            workflow_environment_path("Alignment.yaml")
        container:
            f"{container_base_path}/viroconstrictor_alignment_{get_hash('Alignment')}.sif"
        log:
            f"{logdir}RemoveAdapters_p1_" "{Virus}.{RefID}.{sample}.log",
        benchmark:
            f"{logdir}{bench}RemoveAdapters_p1_" "{Virus}.{RefID}.{sample}.txt"
        threads: config["threads"]["Alignments"]
        resources:
            mem_mb=medium_memory_job,
            runtime=medium_runtime_job,
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
                parameter_name=f"Minimap2_RawAlignParams_{config['platform']}",
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
            minimap2 --cs {params.mm2_alignment_preset} {params.minimap2_base_setting} {params.minimap2_extra_setting} {params.minimap2_alignmentparams} -t {params.mapthreads} {input.ref} {input.fq} 2>> {log} |\
            samtools view -@ {threads} {params.samtools_standard_filters} {params.samtools_extra_filters} -uS 2>> {log} |\
            samtools sort -o {output.bam} >> {log} 2>&1
            samtools index {output.bam} >> {log} 2>&1
            """


if config["platform"] == "illumina" and config["unidirectional"] is False:
    base_mm2_preset = "-ax sr"

    rule remove_adapters_p1:
        input:
            ref=rules.prepare_refs.output,
            fq1=lambda wildcards: SAMPLES[wildcards.sample]["R1"],
            fq2=lambda wildcards: SAMPLES[wildcards.sample]["R2"],
        output:
            bam=f"{datadir}{wc_folder}{cln}{raln}" "{sample}.bam",
            index=f"{datadir}{wc_folder}{cln}{raln}" "{sample}.bam.bai",
        conda:
            workflow_environment_path("Alignment.yaml")
        container:
            f"{container_base_path}/viroconstrictor_alignment_{get_hash('Alignment')}.sif"
        log:
            f"{logdir}" "RemoveAdapters_p1_{Virus}.{RefID}.{sample}.log",
        benchmark:
            f"{logdir}{bench}" "RemoveAdapters_p1_{Virus}.{RefID}.{sample}.log"
        threads: config["threads"]["Alignments"]
        resources:
            mem_mb=medium_memory_job,
            runtime=medium_runtime_job,
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
                parameter_name=f"Minimap2_RawAlignParams_{config['platform']}",
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
            minimap2 --cs {params.mm2_alignment_preset} {params.minimap2_base_setting} {params.minimap2_extra_setting} {params.minimap2_alignmentparams} -t {params.mapthreads} {input.ref} {input.fq1} {input.fq2} 2>> {log} |\
            samtools view -@ {threads} {params.samtools_standard_filters} {params.samtools_extra_filters} -uS 2>> {log} |\
            samtools sort -o {output.bam} >> {log} 2>&1
            samtools index {output.bam} >> {log} 2>&1
            """


rule remove_adapters_p2:
    input:
        rules.remove_adapters_p1.output.bam,
    output:
        f"{datadir}{wc_folder}{cln}{noad}" "{sample}.fastq",
    conda:
        workflow_environment_path("core_scripts.yaml")
    container:
        f"{container_base_path}/viroconstrictor_core_scripts_{get_hash('core_scripts')}.sif"
    log:
        f"{logdir}RemoveAdapters_p2_" "{Virus}.{RefID}.{sample}.log",
    benchmark:
        f"{logdir}{bench}RemoveAdapters_p2_" "{Virus}.{RefID}.{sample}.txt"
    threads: config["threads"]["AdapterRemoval"]
    resources:
        mem_mb=low_memory_job,
        runtime=medium_runtime_job,
    params:
        script="-m main.scripts.clipper",
        pythonpath = f'{Path(workflow.basedir).parent}',
        clipper_filterparams=lambda wc: get_preset_parameter(
            preset_name=SAMPLES[wc.sample]["PRESET"],
            parameter_name=f"Clipper_FilterParams_{config['platform']}",
        ),
    shell:
        """
        PYTHONPATH={params.pythonpath} \
        python {params.script} \
        --input {input} \
        --output {output} \
        {params.clipper_filterparams} \
        --threads {threads}
        """
