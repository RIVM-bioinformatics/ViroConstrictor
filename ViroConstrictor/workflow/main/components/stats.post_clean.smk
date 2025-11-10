
rule qc_clean:
    input:
        rules.qc_filter.output.fq,
    output:
        html=f"{datadir}{wc_folder}{qc_post}" "{sample}_fastqc.html",
        zip=f"{datadir}{wc_folder}{qc_post}" "{sample}_fastqc.zip",
    conda:
        workflow_environment_path("core_scripts.yaml")
    container:
        f"{container_base_path}/viroconstrictor_core_scripts_{get_hash('core_scripts')}.sif"
    log:
        f"{logdir}QC_clean_" "{Virus}.{RefID}.{sample}.log",
    benchmark:
        f"{logdir}{bench}QC_clean_" "{Virus}.{RefID}.{sample}.txt"
    threads: config["threads"]["QC"]
    resources:
        mem_mb=low_memory_job,
        runtime=medium_runtime_job,
    params:
        outdir=f"{datadir}{wc_folder}{qc_post}",
        script=workflow_script_path("scripts/fastqc.sh")
    shell:
        """
        bash {params.script} {input} {params.outdir} {output.html} {output.zip} {log}
        """


def construct_MultiQC_input(_wildcards):
    pre = []
    if config["platform"] in ["nanopore", "iontorrent"]:
        pre = expand(f"{datadir}{qc_pre}" "{sample}_fastqc.zip", sample=SAMPLES)
    elif config["platform"] == "illumina":
        if config["unidirectional"] is True:
            pre.extend(
                expand(f"{datadir}{qc_pre}" "{sample}_fastqc.zip", sample=SAMPLES)
            )
        else:
            pre.extend(
                expand(
                    f"{datadir}{qc_pre}" "{sample}_{read}_fastqc.zip",
                    sample=SAMPLES,
                    read=["R1", "R2"],
                )
            )
    else:
        raise ValueError(
            f"Platform {config['platform']} not recognised. Choose one of [illumina, nanopore, iontorrent]."
        )

    post = expand(
        f"{datadir}{wc_folder}{qc_post}" "{sample}_fastqc.zip",
        zip,
        RefID=p_space.RefID,
        Virus=p_space.Virus,
        sample=p_space.dataframe["sample"],
    )

    fastp_out = expand(
        f"{datadir}{wc_folder}{cln}{qcfilt}{json}" "{sample}_fastp.json",
        zip,
        RefID=p_space.RefID,
        Virus=p_space.Virus,
        sample=p_space.dataframe["sample"],
    )

    return pre + list(post) + fastp_out


rule multiqc_report:
    input:
        construct_MultiQC_input,
    output:
        f"{res}multiqc.html",
        expand(f"{res}{mqc_data}multiqc_" "{program}.txt", program="fastqc"),
    conda:
        workflow_environment_path("Clean.yaml")
    container:
        f"{container_base_path}/viroconstrictor_clean_{get_hash('Clean')}.sif"
    log:
        f"{logdir}MultiQC_report.log",
    benchmark:
        f"{logdir}{bench}MultiQC_report.txt"
    resources:
        mem_mb=high_memory_job,
        runtime=medium_runtime_job,
    params:
        conffile=workflow_script_path("configs/multiqc.yaml"),
        outdir=res,
    shell:
        """
        multiqc -d --force --config {params.conffile} -o {params.outdir} -n multiqc.html {input} > {log} 2>&1
        """
