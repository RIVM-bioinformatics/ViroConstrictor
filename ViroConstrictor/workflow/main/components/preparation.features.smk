
rule prepare_gffs:
    input:
        feats=lambda wc: file if (file := SAMPLES[wc.sample]["FEATURES"]) else "",
        ref=rules.prepare_refs.output,
    output:
        gff=f"{datadir}{wc_folder}{features}" "{sample}_features.gff",
    log:
        f"{logdir}prepare_gffs_" "{Virus}.{RefID}.{sample}.log",
    conda:
        workflow_environment_path("core_scripts.yaml")
    container:
        f"{container_base_path}/viroconstrictor_core_scripts_{get_hash('core_scripts')}.sif"
    resources:
        mem_mb=low_memory_job,
        runtime=low_runtime_job,
    params:
        pythonpath = f'{Path(workflow.basedir).parent}',
        script="-m main.scripts.extract_gff"
    shell:
        """
        PYTHONPATH={params.pythonpath} \
        python {params.script} \
        --input {input.feats} \
        --output {output.gff} \
        --ref_id {wildcards.RefID}
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
    threads: config["threads"]["Index"]
    conda:
        workflow_environment_path("ORF_analysis.yaml")
    container:
        f"{container_base_path}/viroconstrictor_orf_analysis_{get_hash('ORF_analysis')}.sif"
    resources:
        mem_mb=medium_memory_job,
        runtime=medium_runtime_job,
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
