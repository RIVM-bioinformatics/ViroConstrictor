if config["platform"] in ["nanopore", "iontorrent"] or (
    config["platform"] == "illumina" and config["unidirectional"] is True
):
    base_mm2_preset = (
        "-ax sr" if config["platform"] in ["iontorrent", "illumina"] else "-ax map-ont"
    )

    rule qc_raw:
        input:
            lambda wc: SAMPLES[wc.sample]["INPUTFILE"],
        output:
            html=f"{datadir}{qc_pre}" "{sample}_fastqc.html",
            zip=f"{datadir}{qc_pre}" "{sample}_fastqc.zip",
        conda:
            workflow_environment_path("Scripts.yaml")
        container:
            f"{container_base_path}/viroconstrictor_scripts_{get_hash('Scripts')}.sif"
        log:
            f"{logdir}QC_raw_data_" "{sample}.log",
        benchmark:
            f"{logdir}{bench}QC_raw_data_" "{sample}.txt"
        threads: config["threads"]["QC"]
        resources:
            mem_mb=low_memory_job,
            runtime=55,
        params:
            output_dir=f"{datadir}{qc_pre}",
            script=(
                workflow_script_path("scripts/fastqc.sh")
                if (
                    DeploymentMethod.CONDA
                    in workflow.deployment_settings.deployment_method
                )
                is True
                else "/scripts/fastqc.sh"
            ),
        shell:
            """
            bash {params.script} {input} {params.output_dir} {output.html} {output.zip} {log}
            """


if config["platform"] == "illumina" and config["unidirectional"] is False:
    base_mm2_preset = "-ax sr"

    rule qc_raw:
        input:
            lambda wildcards: SAMPLES[wildcards.sample][wildcards.read],
        output:
            html=f"{datadir}{qc_pre}" "{sample}_{read}_fastqc.html",
            zip=f"{datadir}{qc_pre}" "{sample}_{read}_fastqc.zip",
        conda:
            workflow_environment_path("Scripts.yaml")
        container:
            f"{container_base_path}/viroconstrictor_scripts_{get_hash('Scripts')}.sif"
        log:
            f"{logdir}" "QC_raw_data_{sample}_{read}.log",
        benchmark:
            f"{logdir}{bench}" "QC_raw_data_{sample}_{read}.txt"
        threads: config["threads"]["QC"]
        resources:
            mem_mb=lambda: low_memory_job(),
            runtime=55,
        params:
            output_dir=f"{datadir}{qc_pre}",
            script=(
                workflow_script_path("wrappers/fastqc_wrapper.sh")
                if (
                    DeploymentMethod.CONDA
                    in workflow.deployment_settings.deployment_method
                )
                is True
                else "/wrappers/fastqc_wrapper.sh"
            ),
        shell:
            """
            bash {params.script} {input} {params.output_dir} {output.html} {output.zip} {log}
            """
