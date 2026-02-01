# Combined results component for easily accessible aggregate results
# This component creates an additional results structure organized by sample and by virus
# rather than by Virus/RefID, making it easier to navigate when working with
# many references per sample.

# ===========================
# BY SAMPLE AGGREGATION RULES
# ===========================

# Rule to combine consensus sequences by sample with modified headers
rule combine_consensus_by_sample:
    input:
        lambda wc: [
            f"{datadir}Virus~{row['Virus']}/RefID~{row['RefID']}/{cons}{seqs}{wc.sample}.fa"
            for _, row in samples_df[samples_df["sample"] == wc.sample].iterrows()
        ]
    output:
        f"{res}{combined}{by_sample}" "{sample}/consensus.fasta"
    resources:
        mem_mb=low_memory_job,
        runtime=low_runtime_job,
    conda:
        workflow_environment_path("core_scripts.yaml")
    container:
        f"{container_base_path}/viroconstrictor_core_scripts_{get_hash('core_scripts')}.sif"
    log:
        f"{logdir}combine_consensus_by_sample_" "{sample}.log",
    params:
        script="-m main.scripts.combine_fasta",
        pythonpath=f'{Path(workflow.basedir).parent}',
        virus_list=lambda wc: " ".join([
            row['Virus'] for _, row in samples_df[samples_df["sample"] == wc.sample].iterrows()
        ]),
        refid_list=lambda wc: " ".join([
            row['RefID'] for _, row in samples_df[samples_df["sample"] == wc.sample].iterrows()
        ]),
    shell:
        """
        PYTHONPATH={params.pythonpath} \
        python {params.script} \
        --input "empty" \
        --input_files {input} \
        --virus_list {params.virus_list} \
        --refid_list {params.refid_list} \
        --output {output} >> {log} 2>&1
        """


# Rule to combine mutations by sample
rule combine_mutations_by_sample:
    input:
        lambda wc: [
            f"{datadir}Virus~{row['Virus']}/RefID~{row['RefID']}/{aln}{vf}{wc.sample}.tsv"
            for _, row in samples_df[samples_df["sample"] == wc.sample].iterrows()
        ]
    output:
        f"{res}{combined}{by_sample}" "{sample}/mutations.tsv"
    resources:
        mem_mb=low_memory_job,
        runtime=low_runtime_job,
    conda:
        workflow_environment_path("core_scripts.yaml")
    container:
        f"{container_base_path}/viroconstrictor_core_scripts_{get_hash('core_scripts')}.sif"
    log:
        f"{logdir}combine_mutations_by_sample_" "{sample}.log",
    params:
        script="-m main.scripts.combine_tabular",
        pythonpath=f'{Path(workflow.basedir).parent}',
        virus_list=lambda wc: " ".join([
            row['Virus'] for _, row in samples_df[samples_df["sample"] == wc.sample].iterrows()
        ]),
        refid_list=lambda wc: " ".join([
            row['RefID'] for _, row in samples_df[samples_df["sample"] == wc.sample].iterrows()
        ]),
    shell:
        """
        PYTHONPATH={params.pythonpath} \
        python {params.script} \
        --input "empty" \
        --input_files {input} \
        --virus_list {params.virus_list} \
        --refid_list {params.refid_list} \
        --file_type mutations \
        --separator $'\\t' \
        --output {output} >> {log} 2>&1
        """


# Rule to combine width of coverage by sample
rule combine_coverage_by_sample:
    input:
        lambda wc: [
            f"{datadir}Virus~{row['Virus']}/RefID~{row['RefID']}/{boc}{wc.sample}.tsv"
            for _, row in samples_df[samples_df["sample"] == wc.sample].iterrows()
        ]
    output:
        f"{res}{combined}{by_sample}" "{sample}/Width_of_coverage.tsv"
    resources:
        mem_mb=low_memory_job,
        runtime=low_runtime_job,
    conda:
        workflow_environment_path("core_scripts.yaml")
    container:
        f"{container_base_path}/viroconstrictor_core_scripts_{get_hash('core_scripts')}.sif"
    log:
        f"{logdir}combine_coverage_by_sample_" "{sample}.log",
    params:
        script="-m main.scripts.combine_tabular",
        pythonpath=f'{Path(workflow.basedir).parent}',
        virus_list=lambda wc: " ".join([
            row['Virus'] for _, row in samples_df[samples_df["sample"] == wc.sample].iterrows()
        ]),
        refid_list=lambda wc: " ".join([
            row['RefID'] for _, row in samples_df[samples_df["sample"] == wc.sample].iterrows()
        ]),
    shell:
        """
        PYTHONPATH={params.pythonpath} \
        python {params.script} \
        --input "empty" \
        --input_files {input} \
        --virus_list {params.virus_list} \
        --refid_list {params.refid_list} \
        --file_type coverage \
        --separator $'\\t' \
        --output {output} >> {log} 2>&1
        """


# Rule to combine amplicon coverage by sample
rule combine_amplicon_coverage_by_sample:
    input:
        lambda wc: [
            f"{datadir}Virus~{row['Virus']}/RefID~{row['RefID']}/{prim}{wc.sample}_ampliconcoverage.csv"
            for _, row in samples_df[samples_df["sample"] == wc.sample].iterrows()
        ]
    output:
        f"{res}{combined}{by_sample}" "{sample}/Amplicon_coverage.csv"
    resources:
        mem_mb=low_memory_job,
        runtime=low_runtime_job,
    conda:
        workflow_environment_path("core_scripts.yaml")
    container:
        f"{container_base_path}/viroconstrictor_core_scripts_{get_hash('core_scripts')}.sif"
    log:
        f"{logdir}combine_amplicon_coverage_by_sample_" "{sample}.log",
    params:
        script="-m main.scripts.combine_tabular",
        pythonpath=f'{Path(workflow.basedir).parent}',
        virus_list=lambda wc: " ".join([
            row['Virus'] for _, row in samples_df[samples_df["sample"] == wc.sample].iterrows()
        ]),
        refid_list=lambda wc: " ".join([
            row['RefID'] for _, row in samples_df[samples_df["sample"] == wc.sample].iterrows()
        ]),
    shell:
        """
        PYTHONPATH={params.pythonpath} \
        python {params.script} \
        --input "empty" \
        --input_files {input} \
        --virus_list {params.virus_list} \
        --refid_list {params.refid_list} \
        --file_type amplicon_coverage \
        --separator , \
        --output {output} >> {log} 2>&1
        """


# Rule to combine aminoacids by sample when input files exist (one file per feature)
rule combine_aminoacids_by_sample:
    input:
        lambda wc: list(set(
            f"{res}Virus~{row['Virus']}/RefID~{row['RefID']}/{amino}{wc.feature}.faa"
            for _, row in samples_df[samples_df["sample"] == wc.sample].iterrows()
            if (
                pd.notna(row.get("AA_FEAT_NAMES")) and
                isinstance(row["AA_FEAT_NAMES"], (list, tuple)) and
                wc.feature in row["AA_FEAT_NAMES"]
            )
        ))
    output:
        f"{res}{combined}{by_sample}" "{sample}/aminoacids/{feature}.faa"
    resources:
        mem_mb=low_memory_job,
        runtime=low_runtime_job,
    conda:
        workflow_environment_path("core_scripts.yaml")
    container:
        f"{container_base_path}/viroconstrictor_core_scripts_{get_hash('core_scripts')}.sif"
    log:
        f"{logdir}combine_aminoacids_by_sample_" "{sample}.{feature}.log",
    params:
        script="-m main.scripts.extract_sample_from_fasta",
        pythonpath=f'{Path(workflow.basedir).parent}',
    shell:
        """
        mkdir -p $(dirname {output})
        PYTHONPATH={params.pythonpath} \
        python {params.script} \
        --input "empty" \
        --input_files {input} \
        --sample_name {wildcards.sample} \
        --output {output} >> {log} 2>&1
        """


# Rule to create empty aminoacid file when no input files exist
rule combine_aminoacids_by_sample_empty:
    output:
        temp(f"{res}{combined}{by_sample}" "{sample}/aminoacids/{feature}.empty.faa")
    resources:
        mem_mb=low_memory_job,
        runtime=low_runtime_job,
    log:
        f"{logdir}combine_aminoacids_by_sample_empty_" "{sample}.{feature}.log",
    shell:
        """
        mkdir -p $(dirname {output})
        echo "No input files found for sample {wildcards.sample} and feature {wildcards.feature}" >> {log}
        touch {output}
        """


ruleorder: combine_aminoacids_by_sample > combine_aminoacids_by_sample_empty


# ===========================
# BY VIRUS AGGREGATION RULES
# ===========================

# Rule to combine consensus sequences by virus with modified headers
rule combine_consensus_by_virus:
    input:
        lambda wc: [
            f"{datadir}Virus~{wc.Virus}/RefID~{row['RefID']}/{cons}{seqs}{row['sample']}.fa"
            for _, row in samples_df[samples_df["Virus"] == wc.Virus].iterrows()
        ]
    output:
        f"{res}Virus~{{Virus}}/{combined}consensus.fasta"
    resources:
        mem_mb=low_memory_job,
        runtime=low_runtime_job,
    conda:
        workflow_environment_path("core_scripts.yaml")
    container:
        f"{container_base_path}/viroconstrictor_core_scripts_{get_hash('core_scripts')}.sif"
    log:
        f"{logdir}combine_consensus_by_virus_" "{Virus}.log",
    params:
        script="-m main.scripts.combine_fasta",
        pythonpath=f'{Path(workflow.basedir).parent}',
        virus_list=lambda wc: " ".join([
            row['Virus'] for _, row in samples_df[samples_df["Virus"] == wc.Virus].iterrows()
        ]),
        refid_list=lambda wc: " ".join([
            row['RefID'] for _, row in samples_df[samples_df["Virus"] == wc.Virus].iterrows()
        ]),
    shell:
        """
        PYTHONPATH={params.pythonpath} \
        python {params.script} \
        --input "empty" \
        --input_files {input} \
        --virus_list {params.virus_list} \
        --refid_list {params.refid_list} \
        --output {output} >> {log} 2>&1
        """


# Rule to combine mutations by virus
rule combine_mutations_by_virus:
    input:
        lambda wc: [
            f"{datadir}Virus~{wc.Virus}/RefID~{row['RefID']}/{aln}{vf}{row['sample']}.tsv"
            for _, row in samples_df[samples_df["Virus"] == wc.Virus].iterrows()
        ]
    output:
        f"{res}Virus~{{Virus}}/{combined}mutations.tsv"
    resources:
        mem_mb=medium_memory_job,
        runtime=low_runtime_job,
    conda:
        workflow_environment_path("core_scripts.yaml")
    container:
        f"{container_base_path}/viroconstrictor_core_scripts_{get_hash('core_scripts')}.sif"
    log:
        f"{logdir}combine_mutations_by_virus_" "{Virus}.log",
    params:
        script="-m main.scripts.combine_tabular",
        pythonpath=f'{Path(workflow.basedir).parent}',
        virus_list=lambda wc: " ".join([
            row['Virus'] for _, row in samples_df[samples_df["Virus"] == wc.Virus].iterrows()
        ]),
        refid_list=lambda wc: " ".join([
            row['RefID'] for _, row in samples_df[samples_df["Virus"] == wc.Virus].iterrows()
        ]),
    shell:
        """
        PYTHONPATH={params.pythonpath} \
        python {params.script} \
        --input "empty" \
        --input_files {input} \
        --virus_list {params.virus_list} \
        --refid_list {params.refid_list} \
        --file_type mutations \
        --separator $'\\t' \
        --output {output} >> {log} 2>&1
        """


# Rule to combine coverage by virus
rule combine_coverage_by_virus:
    input:
        lambda wc: [
            f"{datadir}Virus~{wc.Virus}/RefID~{row['RefID']}/{boc}{row['sample']}.tsv"
            for _, row in samples_df[samples_df["Virus"] == wc.Virus].iterrows()
        ]
    output:
        f"{res}Virus~{{Virus}}/{combined}Width_of_coverage.tsv"
    resources:
        mem_mb=medium_memory_job,
        runtime=low_runtime_job,
    conda:
        workflow_environment_path("core_scripts.yaml")
    container:
        f"{container_base_path}/viroconstrictor_core_scripts_{get_hash('core_scripts')}.sif"
    log:
        f"{logdir}combine_coverage_by_virus_" "{Virus}.log",
    params:
        script="-m main.scripts.combine_tabular",
        pythonpath=f'{Path(workflow.basedir).parent}',
        virus_list=lambda wc: " ".join([
            row['Virus'] for _, row in samples_df[samples_df["Virus"] == wc.Virus].iterrows()
        ]),
        refid_list=lambda wc: " ".join([
            row['RefID'] for _, row in samples_df[samples_df["Virus"] == wc.Virus].iterrows()
        ]),
    shell:
        """
        PYTHONPATH={params.pythonpath} \
        python {params.script} \
        --input "empty" \
        --input_files {input} \
        --virus_list {params.virus_list} \
        --refid_list {params.refid_list} \
        --file_type coverage \
        --separator $'\\t' \
        --output {output} >> {log} 2>&1
        """


# Rule to combine amplicon coverage by virus
rule combine_amplicon_coverage_by_virus:
    input:
        lambda wc: [
            f"{datadir}Virus~{wc.Virus}/RefID~{row['RefID']}/{prim}{row['sample']}_ampliconcoverage.csv"
            for _, row in samples_df[samples_df["Virus"] == wc.Virus].iterrows()
        ]
    output:
        f"{res}Virus~{{Virus}}/{combined}Amplicon_coverage.csv"
    resources:
        mem_mb=medium_memory_job,
        runtime=low_runtime_job,
    conda:
        workflow_environment_path("core_scripts.yaml")
    container:
        f"{container_base_path}/viroconstrictor_core_scripts_{get_hash('core_scripts')}.sif"
    log:
        f"{logdir}combine_amplicon_coverage_by_virus_" "{Virus}.log",
    params:
        script="-m main.scripts.combine_tabular",
        pythonpath=f'{Path(workflow.basedir).parent}',
        virus_list=lambda wc: " ".join([
            row['Virus'] for _, row in samples_df[samples_df["Virus"] == wc.Virus].iterrows()
        ]),
        refid_list=lambda wc: " ".join([
            row['RefID'] for _, row in samples_df[samples_df["Virus"] == wc.Virus].iterrows()
        ]),
    shell:
        """
        PYTHONPATH={params.pythonpath} \
        python {params.script} \
        --input "empty" \
        --input_files {input} \
        --virus_list {params.virus_list} \
        --refid_list {params.refid_list} \
        --file_type amplicon_coverage \
        --separator , \
        --output {output} >> {log} 2>&1
        """


# Rule to combine aminoacids by virus (one file per feature)
rule combine_aminoacids_by_virus:
    input:
        lambda wc: list(set([
            f"{res}Virus~{wc.Virus}/RefID~{row['RefID']}/{amino}{wc.feature}.faa"
            for _, row in samples_df[samples_df["Virus"] == wc.Virus].iterrows()
            if (
                pd.notna(row.get("AA_FEAT_NAMES")) and
                isinstance(row["AA_FEAT_NAMES"], (list, tuple)) and 
                wc.feature in row["AA_FEAT_NAMES"]
            )
        ]))
    output:
        f"{res}Virus~{{Virus}}/{combined}aminoacids/{{feature}}.faa"
    resources:
        mem_mb=low_memory_job,
        runtime=low_runtime_job,
    log:
        f"{logdir}combine_aminoacids_by_virus_" "{Virus}.{feature}.log",
    shell:
        """
        mkdir -p $(dirname {output})
        cat {input} > {output} 2>> {log}
        """


# Rule to create empty aminoacid file when no input files exist for virus
rule combine_aminoacids_by_virus_empty:
    output:
        temp(f"{res}Virus~{{Virus}}/{combined}aminoacids/{{feature}}.empty.faa")
    resources:
        mem_mb=low_memory_job,
        runtime=low_runtime_job,
    log:
        f"{logdir}combine_aminoacids_by_virus_empty_" "{Virus}.{feature}.log",
    shell:
        """
        mkdir -p $(dirname {output:q})
        echo "No input files found for virus {wildcards.Virus} and feature {wildcards.feature}" >> {log:q}
        touch {output:q}
        """


ruleorder: combine_aminoacids_by_virus > combine_aminoacids_by_virus_empty


# =============================
# ALL SAMPLES AGGREGATION RULES
# =============================

# Global aggregation rule for consensus (from by_sample results)
rule combine_all_consensus:
    input:
        expand(
            f"{res}{combined}{by_sample}" "{sample}/consensus.fasta",
            sample=samples_df["sample"].unique()
        )
    output:
        f"{res}{combined}{all_samples}all_consensus.fasta"
    resources:
        mem_mb=low_memory_job,
        runtime=low_runtime_job,
    shell:
        """
        cat {input} > {output}
        """


# Global aggregation rule for mutations (from by_sample results)
rule combine_all_mutations:
    input:
        expand(
            f"{res}{combined}{by_sample}" "{sample}/mutations.tsv",
            sample=samples_df["sample"].unique()
        )
    output:
        f"{res}{combined}{all_samples}all_mutations.tsv"
    resources:
        mem_mb=medium_memory_job,
        runtime=low_runtime_job,
    conda:
        workflow_environment_path("core_scripts.yaml")
    container:
        f"{container_base_path}/viroconstrictor_core_scripts_{get_hash('core_scripts')}.sif"
    log:
        f"{logdir}combine_all_mutations.log",
    params:
        script="-m main.scripts.aggregate_combined_files",
        pythonpath=f'{Path(workflow.basedir).parent}',
    shell:
        """
        PYTHONPATH={params.pythonpath} \
        python {params.script} \
        --input "empty" \
        --input_files {input} \
        --file_type mutations \
        --separator $'\\t' \
        --output {output} >> {log} 2>&1
        """


# Global aggregation rule for coverage (from by_sample results)
rule combine_all_coverage:
    input:
        expand(
            f"{res}{combined}{by_sample}" "{sample}/Width_of_coverage.tsv",
            sample=samples_df["sample"].unique()
        )
    output:
        f"{res}{combined}{all_samples}all_width_of_coverage.tsv"
    resources:
        mem_mb=medium_memory_job,
        runtime=low_runtime_job,
    conda:
        workflow_environment_path("core_scripts.yaml")
    container:
        f"{container_base_path}/viroconstrictor_core_scripts_{get_hash('core_scripts')}.sif"
    log:
        f"{logdir}combine_all_coverage.log",
    params:
        script="-m main.scripts.aggregate_combined_files",
        pythonpath=f'{Path(workflow.basedir).parent}',
    shell:
        """
        PYTHONPATH={params.pythonpath} \
        python {params.script} \
        --input "empty" \
        --input_files {input:q} \
        --file_type coverage \
        --separator $'\\t' \
        --output {output} >> {log} 2>&1
        """


# Global aggregation rule for amplicon coverage (from by_sample results)
rule combine_all_amplicon_coverage:
    input:
        expand(
            f"{res}{combined}{by_sample}" "{sample}/Amplicon_coverage.csv",
            sample=samples_df["sample"].unique()
        )
    output:
        f"{res}{combined}{all_samples}all_amplicon_coverage.csv"
    resources:
        mem_mb=medium_memory_job,
        runtime=low_runtime_job,
    conda:
        workflow_environment_path("core_scripts.yaml")
    container:
        f"{container_base_path}/viroconstrictor_core_scripts_{get_hash('core_scripts')}.sif"
    log:
        f"{logdir}combine_all_amplicon_coverage.log",
    params:
        script="-m main.scripts.aggregate_combined_files",
        pythonpath=f'{Path(workflow.basedir).parent}',
    shell:
        """
        PYTHONPATH={params.pythonpath} \
        python {params.script} \
        --input "empty" \
        --input_files {input:q} \
        --file_type amplicon_coverage \
        --separator , \
        --output {output} >> {log} 2>&1
        """

# Global aggregation rule for aminoacids (from by_sample results, one file per feature)
rule combine_all_aminoacids:
    input:
        lambda wc: [
            f"{res}{combined}{by_sample}{sample}/aminoacids/{wc.feature}.faa"
            for sample in samples_df["sample"].unique()
            if any(
                pd.notna(row.get("AA_FEAT_NAMES")) and
                isinstance(row["AA_FEAT_NAMES"], (list, tuple)) and
                wc.feature in row["AA_FEAT_NAMES"]
                for _, row in samples_df[samples_df["sample"] == sample].iterrows()
            )
        ]
    output:
        f"{res}{combined}{all_samples}aminoacids/{{feature}}.faa"
    resources:
        mem_mb=low_memory_job,
        runtime=low_runtime_job,
    log:
        f"{logdir}combine_all_aminoacids_" "{feature}.log",
    shell:
        """
        mkdir -p $(dirname {output})
        cat {input} > {output} 2>> {log}
        """


# Rule to create empty aminoacid file when no input files exist
rule combine_all_aminoacids_empty:
    output:
        temp(f"{res}{combined}{all_samples}aminoacids/{{feature}}.empty.faa")
    resources:
        mem_mb=low_memory_job,
        runtime=low_runtime_job,
    log:
        f"{logdir}combine_all_aminoacids_empty_" "{feature}.log",
    shell:
        """
        mkdir -p $(dirname {output:q})
        echo "No input files found for feature {wildcards.feature}" >> {log:q}
        touch {output:q}
        """


ruleorder: combine_all_aminoacids > combine_all_aminoacids_empty