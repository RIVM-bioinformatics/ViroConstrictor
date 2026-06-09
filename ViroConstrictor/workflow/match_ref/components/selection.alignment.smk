if config["unidirectional"] is True:
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
            {params.align_bin} {params.align_flags} -t {params.mapthreads} {input.ref} {input.fq} 2>> {log} |\
            {params.filter_bin} view -@ {threads} {params.filter_flags} -uS 2>> {log} |\
            {params.filter_bin} sort -o {output.bam} >> {log} 2>&1
            {params.filter_bin} index {output.bam} >> {log} 2>&1
            """


else:
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
        runtime=low_runtime_job,
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
