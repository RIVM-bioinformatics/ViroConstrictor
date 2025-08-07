
rule filter_fasta2bed:
    input:
        ref=rules.group_and_rename_refs.output.groupedrefs,
        refdata=rules.filter_gff.output.groupedstats,
        prm=lambda wc: (
            ""
            if SAMPLES[wc.sample]["PRIMERS"].endswith(".bed")
            else SAMPLES[wc.sample]["PRIMERS"]
        ),
    output:
        bed=f"{datadir}{matchref}" "{sample}_primers.bed",
        groupedstats=temp(f"{datadir}{matchref}" "{sample}_step2.csv"),
    conda:
        workflow_environment_path("Clean.yaml")
    container:
        f"{container_base_path}/viroconstrictor_clean_{get_hash('Clean')}.sif"
    threads: 1
    resources:
        mem_mb=low_memory_job,
        runtime=55,
    log:
        f"{logdir}Fasta2Bed_" "{sample}.log",
    params:
        pr_mm_rate=lambda wc: SAMPLES[wc.sample]["PRIMER-MISMATCH-RATE"],
    shell:
        """
        python -m AmpliGone.fasta2bed \
            --primers {input.prm} \
            --reference {input.ref} \
            --output {output.bed} \
            --primer-mismatch-rate {params.pr_mm_rate} > {log}
        awk -F ',' -v OFS=',' '{{ $(NF+1) = (NR==1 ? "Primer_file" : "{output.bed}"); print }}' {input.refdata} > {output.groupedstats}
        """


ruleorder: filter_fasta2bed > filter_bed > touch_primers


rule filter_bed:
    input:
        prm=lambda wc: SAMPLES[wc.sample]["PRIMERS"],
        refdata=rules.filter_gff.output.groupedstats,
    output:
        bed=f"{datadir}{matchref}" "{sample}_primers.bed",
        groupedstats=temp(f"{datadir}{matchref}" "{sample}_step2.csv"),
    threads: 1
    resources:
        mem_mb=low_memory_job,
        runtime=55,
    log:
        f"{logdir}FilterBed_" "{sample}.log",
    params:
        script=(
            workflow_script_path("script/filter_bed.py")
            if (
                DeploymentMethod.CONDA
                in workflow.deployment_settings.deployment_method
            )
            is True
            else "/match_ref_scripts/filter_bed.py"
        ),
    conda:
        workflow_environment_path("Scripts.yaml")
    container:
        f"{container_base_path}/viroconstrictor_scripts_{get_hash('Scripts')}.sif"
    shell:
        """
        python {params.script} {input.prm} {input.refdata} {output.bed} {output.groupedstats}
        """


rule touch_primers:
    input:
        refdata=rules.filter_gff.output.groupedstats,
    output:
        bed=f"{datadir}{matchref}" "{sample}_primers.bed",
        groupedstats=temp(f"{datadir}{matchref}" "{sample}_step2.csv"),
    threads: 1
    resources:
        mem_mb=low_memory_job,
        runtime=55,
    shell:
        """
        cp {input.refdata} {output.groupedstats}
        touch {output.bed}
        """
