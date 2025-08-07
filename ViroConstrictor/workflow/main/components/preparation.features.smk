
rule prepare_gffs:
    input:
        feats=lambda wc: file if (file := SAMPLES[wc.sample]["FEATURES"]) else "",
        ref=rules.prepare_refs.output,
    output:
        gff=f"{datadir}{wc_folder}{features}" "{sample}_features.gff",
    log:
        f"{logdir}prepare_gffs_" "{Virus}.{RefID}.{sample}.log",
    benchmark:
        f"{logdir}{bench}prepare_gffs_" "{Virus}.{RefID}.{sample}.txt"
    conda:
        workflow_environment_path("Scripts.yaml")
    container:
        f"{container_base_path}/viroconstrictor_scripts_{get_hash('Scripts')}.sif"
    resources:
        mem_mb=low_memory_job,
        runtime=55,
    params:
        script=(
            workflow_script_path("scripts/extract_gff.py")
            if (
                DeploymentMethod.CONDA
                in workflow.deployment_settings.deployment_method
            )
            is True
            else "/scripts/extract_gff.py"
        ),
    shell:
        """
        python {params.script} {input.feats} {output.gff} {wildcards.RefID}
        """


ruleorder: prepare_gffs > prodigal


rule prodigal:
    input:
        ref=rules.prepare_refs.output,
    output:
        gff=f"{datadir}{wc_folder}{features}" "{sample}_features.gff",
        aa=f"{datadir}{wc_folder}{features}" "{sample}_features.aa.fasta",
        nt=f"{datadir}{wc_folder}{features}" "{sample}_features.nt.fasta",
    log:
        f"{logdir}prepare_gffs_" "{Virus}.{RefID}.{sample}.log",
    benchmark:
        f"{logdir}{bench}prepare_gffs_" "{Virus}.{RefID}.{sample}.txt"
    threads: config["threads"]["Index"]
    conda:
        workflow_environment_path("ORF_analysis.yaml")
    container:
        f"{container_base_path}/viroconstrictor_orf_analysis_{get_hash('ORF_analysis')}.sif"
    resources:
        mem_mb=medium_memory_job,
        runtime=55,
    params:
        prodigal_method="meta",
        prodigal_outformat="gff",
    shell:
        """
        prodigal -q -i {input.ref} \
            -a {output.aa} \
            -d {output.nt} \
            -o {output.gff} \
            -p {params.prodigal_method} \
            -f {params.prodigal_outformat} > {log} 2>&1
        """
