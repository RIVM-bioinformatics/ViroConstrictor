rule qc_filter:
    input:
        rules.remove_adapters_p2.output,
    output:
        fq=f"{datadir}{wc_folder}{cln}{qcfilt}" "{sample}.fastq",
        html=f"{datadir}{wc_folder}{cln}{qcfilt}{html}" "{sample}_fastplong.html",
        json=f"{datadir}{wc_folder}{cln}{qcfilt}{json}" "{sample}_fastplong.json",
    conda:
        workflow_environment_path("Clean.yaml")
    container:
        f"{container_base_path}/viroconstrictor_clean_{get_hash('Clean')}.sif"
    log:
        f"{logdir}QC_filter_" "{Virus}.{RefID}.{sample}.log",
    threads: config["threads"]["QC"]
    resources:
        mem_mb=low_memory_job,
        runtime=low_runtime_job,
    params:
        binary = lambda wc: get_binary(
            stage_identifier=VC_STAGE,
            preset_name=SAMPLES[wc.sample]["PRESET"],
            task="qc",
            platform=config["platform"],
        ),
        flags = lambda wc: get_flags(
            stage_identifier=VC_STAGE,
            preset_name=SAMPLES[wc.sample]["PRESET"],
            task="qc",
            platform=config["platform"],
        ),
    shell:
        """
        {params.binary} --thread {threads} -i {input} \
            {params.flags} \
            -o {output.fq} -h {output.html} -j {output.json} > {log} 2>&1
        """