
rule qc_filter:
    input:
        rules.remove_adapters_p2.output,
    output:
        fq=f"{datadir}{wc_folder}{cln}{qcfilt}" "{sample}.fastq",
        html=f"{datadir}{wc_folder}{cln}{qcfilt}{html}" "{sample}_fastp.html",
        json=f"{datadir}{wc_folder}{cln}{qcfilt}{json}" "{sample}_fastp.json",
    conda:
        workflow_environment_path("Clean.yaml")
    container:
        f"{container_base_path}/viroconstrictor_clean_{get_hash('Clean')}.sif"
    log:
        f"{logdir}QC_filter_" "{Virus}.{RefID}.{sample}.log",
    benchmark:
        f"{logdir}{bench}QC_filter_" "{Virus}.{RefID}.{sample}.txt"
    threads: config["threads"]["QC"]
    resources:
        mem_mb=low_memory_job,
        runtime=low_runtime_job,
    params:
        adapter_removal_settings=lambda wc: get_preset_parameter(
            preset_name=SAMPLES[wc.sample]["PRESET"],
            parameter_name=f"FastP_AdapterRemoval_Settings_{config['platform']}",
        ),
        per_read_qc_settings=lambda wc: get_preset_parameter(
            preset_name=SAMPLES[wc.sample]["PRESET"],
            parameter_name=f"Fastp_PerReadQC_Settings_{config['platform']}",
        ),
        slidingwindow_qc_settings=lambda wc: get_preset_parameter(
            preset_name=SAMPLES[wc.sample]["PRESET"],
            parameter_name=f"Fastp_QC_SlidingWindow_Settings_{config['platform']}",
        ),
        lowcomplexity_filter_settings=lambda wc: get_preset_parameter(
            preset_name=SAMPLES[wc.sample]["PRESET"],
            parameter_name=f"Fastp_LowComplexityFilter_Settings_{config['platform']}",
        ),
        minlength_filter_settings=lambda wc: get_preset_parameter(
            preset_name=SAMPLES[wc.sample]["PRESET"],
            parameter_name=f"Fastp_MinReadLength_Settings_{config['platform']}",
        ),
    shell:
        """
        fastp --thread {threads} -i {input} \
            {params.adapter_removal_settings} \
            {params.per_read_qc_settings} \
            {params.slidingwindow_qc_settings} \
            {params.lowcomplexity_filter_settings} \
            {params.minlength_filter_settings} \
            -o {output.fq} -h {output.html} -j {output.json} > {log} 2>&1
        """
