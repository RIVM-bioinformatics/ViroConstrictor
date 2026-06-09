if config["unidirectional"] is True:
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
        threads: config["threads"]["Alignments"]
        resources:
            mem_mb=medium_memory_job,
            runtime=medium_runtime_job,
        params:
            mapthreads=config["threads"]["Alignments"] - 1,
            align_bin = lambda wc: get_binary(
                stage_identifier=VC_STAGE,
                preset_name=SAMPLES[wc.sample]["PRESET"],
                task="raw_alignment",
                platform=config["platform"],
            ),
            align_flags = lambda wc: get_flags(
                stage_identifier=VC_STAGE,
                preset_name=SAMPLES[wc.sample]["PRESET"],
                task="raw_alignment",
                platform=config["platform"],
            ),
            filter_bin = lambda wc: get_binary(
                stage_identifier=VC_STAGE,
                preset_name=SAMPLES[wc.sample]["PRESET"],
                task="alignment_filtering",
                platform=config["platform"],
            ),
            filter_flags = lambda wc: get_flags(
                stage_identifier=VC_STAGE,
                preset_name=SAMPLES[wc.sample]["PRESET"],
                task="alignment_filtering",
                platform=config["platform"],
            ),
        shell:
            """
            {params.align_bin} {params.align_flags} -t {params.mapthreads} {input.ref} {input.fq} 2>> {log} |\
            {params.filter_bin} view -@ {threads} {params.filter_flags} -uS 2>> {log} |\
            {params.filter_bin} sort -o {output.bam} >> {log} 2>&1
            {params.filter_bin} index {output.bam} >> {log} 2>&1
            """


else:
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
        threads: config["threads"]["Alignments"]
        resources:
            mem_mb=medium_memory_job,
            runtime=medium_runtime_job,
        params:
            mapthreads=config["threads"]["Alignments"] - 1,
            align_bin = lambda wc: get_binary(
                stage_identifier=VC_STAGE,
                preset_name=SAMPLES[wc.sample]["PRESET"],
                task="alignment",
                platform=config["platform"],
            ),
            align_flags = lambda wc: get_flags(
                stage_identifier=VC_STAGE,
                preset_name=SAMPLES[wc.sample]["PRESET"],
                task="alignment",
                platform=config["platform"],
            ),
            filter_bin = lambda wc: get_binary(
                stage_identifier=VC_STAGE,
                preset_name=SAMPLES[wc.sample]["PRESET"],
                task="alignment_filtering",
                platform=config["platform"],
            ),
            filter_flags = lambda wc: get_flags(
                stage_identifier=VC_STAGE,
                preset_name=SAMPLES[wc.sample]["PRESET"],
                task="alignment_filtering",
                platform=config["platform"],
            ),
        shell:
            """
            {params.align_bin} {params.align_flags} -t {params.mapthreads} {input.ref} {input.fq1} {input.fq2} 2>> {log} |\
            {params.filter_bin} view -@ {threads} {params.filter_flags} -uS 2>> {log} |\
            {params.filter_bin} sort -o {output.bam} >> {log} 2>&1
            {params.filter_bin} index {output.bam} >> {log} 2>&1
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
    threads: config["threads"]["AdapterRemoval"]
    resources:
        mem_mb=low_memory_job,
        runtime=medium_runtime_job,
    params:
        script="-m main.scripts.clipper",
        pythonpath = f'{Path(workflow.basedir).parent}',
        flags = lambda wc: get_flags(
            stage_identifier=VC_STAGE,
            preset_name=SAMPLES[wc.sample]["PRESET"],
            task="adapter_removal",
            platform=config["platform"],
        ),
    shell:
        """
        PYTHONPATH={params.pythonpath} \
        python {params.script} \
        --input {input} \
        --output {output} \
        {params.flags} \
        --threads {threads} >> {log} 2>&1
        """
