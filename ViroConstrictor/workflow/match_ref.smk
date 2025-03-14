import os
import yaml
import sys
import logging

import pandas as pd

from presets import get_preset_parameter
from containers import get_hash

from Bio import SeqIO, SeqRecord
from snakemake.utils import Paramspace, min_version
import snakemake

from directories import *

from rich import print
min_version("7.15")



logger = logging.getLogger()
logger.handlers.clear()

from snakemake import logger

# Elevate the log level of all output generated by the snakemake.logging module to CRITICAL in order to suppress it when snakemake is calling itself in a downstream process.
if "--snakefile" in sys.argv:
    logging.getLogger("snakemake.logging").setLevel(logging.CRITICAL)
import ViroConstrictor


SAMPLES = {}
with open(config["sample_sheet"]) as sample_sheet_file:
    SAMPLES = yaml.safe_load(sample_sheet_file)


def read_fasta(fasta_file: str) -> list[SeqRecord.SeqRecord]:
    """
    Read a FASTA file and return a list of SeqRecord objects.

    Parameters
    ----------
    fasta_file : str
        The path to the FASTA file to be read.

    Returns
    -------
    list[SeqRecord.SeqRecord]
        A list of SeqRecord objects representing the sequences in the FASTA file.

    """
    return list(SeqIO.parse(fasta_file, "fasta"))


def segmented_ref_groups(df: pd.DataFrame) -> pd.DataFrame:
    """
    Filter required ref groups from main reference file.

    Parameters
    ----------
    df : pd.DataFrame
        A pandas DataFrame containing the reference file information.

    Returns
    -------
    pd.DataFrame
        A pandas DataFrame containing the filtered reference file information.

    """
    for index, row in df.iterrows():
        # if the value in the "SEGMENTED" column is False then place a None string in the segment column. 
        # Ensure that the value is a string and not a NoneType.
        if not row["SEGMENTED"]:
            df.at[index, "segment"] = {"None"}
            continue
        refs = read_fasta(row["REFERENCE"])
        unique_groups = {x.description.split(" ")[1].split("|")[0] for x in refs}
        if len(unique_groups) < 2:
            df.drop(index, inplace=True)
            continue
        df.at[index, "segment"] = [unique_groups]
    # ensure the segment column is a flat set of strings and not a list with a set inside.
    df["segment"] = df["segment"].apply(lambda x: x.pop() if isinstance(x, list) else x)
    # TODO: Check this later to see if this can be done in the first pass (when the values are being set in the first place) instead of fixing afterwards. Not entirely sure why the values are being generated as a list with a set inside one time and just as a normal set the other time.
    return df


samples_df = (
    pd.DataFrame(SAMPLES)
    .transpose()
    .reset_index()
    .rename(columns=dict(index="sample", VIRUS="Virus"))
)
samples_df = segmented_ref_groups(samples_df)
samples_df = samples_df.explode("segment")
p_space = Paramspace(
    samples_df[["Virus", "segment", "sample"]], filename_params=["sample"]
)
wc_folder = "/".join(p_space.wildcard_pattern.split("/")[:-1]) + "/"




def low_memory_job(wildcards, threads, attempt):
    if config["computing_execution"] == "local":
        return min(attempt * threads * 1 * 1000, config["max_local_mem"])
    return attempt * threads * 1 * 1000


def medium_memory_job(wildcards, threads, attempt):
    if config["computing_execution"] == "local":
        return min(attempt * threads * 2 * 1000, config["max_local_mem"])
    return attempt * threads * 2 * 1000


def high_memory_job(wildcards, threads, attempt):
    if config["computing_execution"] == "local":
        return min(attempt * threads * 4 * 1000, config["max_local_mem"])
    return attempt * threads * 4 * 1000


wildcard_constraints:
    # regular expression to match only alphanumeric characters, underscores, dashes, and dots. exclude '/' and only match the first part of the string.
    segment="[\w\-\.\d]+",
    Virus="[\w\-\.\d]+",
    sample="[\w\-\.\d]+",


localrules:
    all,
    filter_references,
    concat_frames,
    group_and_rename_refs,
    count_mapped_reads,
    filter_bed,
    filter_gff,


rule all:
    input:
        f"{datadir}" "match_ref_results.pkl",
        expand(
            f"{datadir}{matchref}" "{sample}_primers.bed",
            sample=p_space.dataframe["sample"],
        ),
        expand(
            f"{datadir}{matchref}" "{sample}_refs.fasta",
            sample=p_space.dataframe["sample"],
        ),
        expand(
            f"{datadir}{matchref}" "{sample}_feats.gff",
            sample=p_space.dataframe["sample"],
        ),


rule filter_references:
    input:
        lambda wc: SAMPLES[wc.sample]["REFERENCE"],
    output:
        temp(f"{datadir}{matchref}{wc_folder}" "{sample}_refs.fasta"),
    resources:
        mem=low_memory_job,
    threads: 1
    log:
        f"{logdir}prepare_refs" "{Virus}.{segment}.{sample}.log",
    benchmark:
        f"{logdir}{bench}MR_prepare_refs" "{Virus}.{segment}.{sample}.txt"
    conda:
        f"{conda_envs}Scripts.yaml"
    container:
        f"{config['container_cache']}/viroconstrictor_scripts_{get_hash('Scripts')}.sif"
    params:
        script=srcdir("scripts/match_ref/filter_references.py") if config["use-conda"] is True and config["use-singularity"] is False else "/scripts/match_ref/filter_references.py",
    shell:
        """
        python {params.script} {input} {output} {wildcards.segment} >> {log} 2>&1
        """


if config["platform"] in ["nanopore", "iontorrent"]:
    base_mapping_settings = (
        "-ax sr" if config["platform"] == "iontorrent" else "-ax map-ont"
    )

    rule align_to_refs:
        input:
            ref=rules.filter_references.output,
            fq=lambda wc: SAMPLES[wc.sample]["INPUTFILE"],
        output:
            bam=temp(f"{datadir}{matchref}{wc_folder}" "{sample}.bam"),
            index=temp(f"{datadir}{matchref}{wc_folder}" "{sample}.bam.bai"),
        conda:
            f"{conda_envs}Alignment.yaml"
        container:
            f"{config['container_cache']}/viroconstrictor_alignment_{get_hash('Alignment')}.sif"
        log:
            f"{logdir}AlignMR_" "{Virus}.{segment}.{sample}.log",
        benchmark:
            f"{logdir}{bench}AlignMR_" "{Virus}.{segment}.{sample}.txt"
        threads: config["threads"]["Alignments"]
        resources:
            mem=medium_memory_job,

        params:
            mapthreads=config["threads"]["Alignments"] - 1,
            mapping_base_settings=base_mapping_settings,
            mapping_additionalsettings=lambda wc: get_preset_parameter(
                preset_name=SAMPLES[wc.sample]["PRESET"],
                parameter_name="RawAlign_AdditionalSettings",
            ),
            MappingMatchRefExtraSetting=lambda wc: get_preset_parameter(
                preset_name=SAMPLES[wc.sample]["PRESET"],
                parameter_name="MatchRef_AdditionalAlignmentSettings",
            ),
            BAMfilters=lambda wc: get_preset_parameter(
                preset_name=SAMPLES[wc.sample]["PRESET"],
                parameter_name="BaseBAMFilters",
            ),
            AdditionalBamFilter=lambda wc: get_preset_parameter(
                preset_name=SAMPLES[wc.sample]["PRESET"],
                parameter_name="MatchRef_AdditionalBamFilters",
            ),
        shell:
            """
            minimap2 {params.mapping_base_settings} {params.mapping_additionalsettings} {params.MappingMatchRefExtraSetting} -t {params.mapthreads} {input.ref} {input.fq} 2>> {log} |\
            samtools view -@ {threads} {params.BAMfilters} {params.AdditionalBamFilter} -uS 2>> {log} |\
            samtools sort -o {output.bam} >> {log} 2>&1
            samtools index {output.bam} >> {log} 2>&1
            """


if config["platform"] == "illumina":

    rule align_to_refs:
        input:
            ref=rules.filter_references.output,
            fq1=lambda wildcards: SAMPLES[wildcards.sample]["R1"],
            fq2=lambda wildcards: SAMPLES[wildcards.sample]["R2"],
        output:
            bam=temp(f"{datadir}{matchref}{wc_folder}" "{sample}.bam"),
            index=temp(f"{datadir}{matchref}{wc_folder}" "{sample}.bam.bai"),
        conda:
            f"{conda_envs}Alignment.yaml"
        container:
            f"{config['container_cache']}/viroconstrictor_alignment_{get_hash('Alignment')}.sif"
        log:
            f"{logdir}" "AlignMR_{Virus}.{segment}.{sample}.log",
        benchmark:
            f"{logdir}{bench}" "AlignMR_{Virus}.{segment}.{sample}.log"
        threads: config["threads"]["Alignments"]
        resources:
            mem=medium_memory_job,

        params:
            mapthreads=config["threads"]["Alignments"] - 1,
            mapping_additionalsettings=lambda wc: get_preset_parameter(
                preset_name=SAMPLES[wc.sample]["PRESET"],
                parameter_name="RawAlign_AdditionalSettings",
            ),
            MappingMatchRefExtraSetting=lambda wc: get_preset_parameter(
                preset_name=SAMPLES[wc.sample]["PRESET"],
                parameter_name="MatchRef_AdditionalAlignmentSettings",
            ),
            BAMfilters=lambda wc: get_preset_parameter(
                preset_name=SAMPLES[wc.sample]["PRESET"],
                parameter_name="BaseBAMFilters",
            ),
            AdditionalBamFilter=lambda wc: get_preset_parameter(
                preset_name=SAMPLES[wc.sample]["PRESET"],
                parameter_name="MatchRef_AdditionalBamFilters",
            ),
        shell:
            """
            minimap2 -ax sr {params.mapping_additionalsettings} {params.MappingMatchRefExtraSetting} -t {params.mapthreads} {input.ref} {input.fq1} {input.fq2} 2>> {log} |\
            samtools view -@ {threads} {params.BAMfilters} {params.AdditionalBamFilter} -uS 2>> {log} |\
            samtools sort -o {output.bam} >> {log} 2>&1
            samtools index {output.bam} >> {log} 2>&1
            """


rule count_mapped_reads:
    input:
        bam=rules.align_to_refs.output.bam,
        index=rules.align_to_refs.output.index,
    output:
        temp(f"{datadir}{matchref}{wc_folder}" "{sample}_count.csv"),
    conda:
        f"{conda_envs}Scripts.yaml"
    container:
        f"{config['container_cache']}/viroconstrictor_scripts_{get_hash('Scripts')}.sif"
    threads: 1
    resources:
        mem=low_memory_job,
    log:
        f"{logdir}CountMR_" "{Virus}.{segment}.{sample}.log",
    params:
        script=srcdir("scripts/match_ref/count_mapped_reads.py") if config["use-conda"] is True and config["use-singularity"] is False else "/scripts/match_ref/count_mapped_reads.py",
    shell:
        """
        python {params.script} {input.bam} {output} >> {log} 2>&1
        """


rule filter_best_matching_ref:
    input:
        stats=rules.count_mapped_reads.output,
        ref=rules.filter_references.output,
    output:
        filtref=temp(f"{datadir}{matchref}{wc_folder}" "{sample}_best_ref.fasta"),
        filtcount=temp(f"{datadir}{matchref}{wc_folder}" "{sample}_best_ref.csv"),
    conda:
        f"{conda_envs}Scripts.yaml"
    container:
        f"{config['container_cache']}/viroconstrictor_scripts_{get_hash('Scripts')}.sif"
    threads: 1
    resources:
        mem=low_memory_job,
    log:
        f"{logdir}FilterBR_" "{Virus}.{segment}.{sample}.log",
    params:
        script=srcdir("scripts/match_ref/filter_best_matching_ref.py") if config["use-conda"] is True and config["use-singularity"] is False else "/scripts/match_ref/filter_best_matching_ref.py",
    shell:
        """
        python {params.script} {input.stats} {input.ref} {output.filtref} {output.filtcount} >> {log} 2>&1
        """


rule group_and_rename_refs:
    input:
        ref=lambda wildcards: expand(
            f"{datadir}{matchref}{wc_folder}" "{sample}_best_ref.fasta",
            zip,
            Virus=p_space.dataframe.loc[
                p_space.dataframe["sample"] == wildcards.sample, "Virus"
            ],
            segment=p_space.dataframe.loc[
                p_space.dataframe["sample"] == wildcards.sample, "segment"
            ],
            allow_missing=True,
        ),
        stats=lambda wildcards: expand(
            f"{datadir}{matchref}{wc_folder}" "{sample}_best_ref.csv",
            zip,
            Virus=p_space.dataframe.loc[
                p_space.dataframe["sample"] == wildcards.sample, "Virus"
            ],
            segment=p_space.dataframe.loc[
                p_space.dataframe["sample"] == wildcards.sample, "segment"
            ],
            allow_missing=True,
        ),
    output:
        groupedrefs=f"{datadir}{matchref}" "{sample}_refs.fasta",
        groupedstats=temp(f"{datadir}{matchref}" "{sample}_refs.csv"),
    conda:
        f"{conda_envs}Scripts.yaml"
    container:
        f"{config['container_cache']}/viroconstrictor_scripts_{get_hash('Scripts')}.sif"
    threads: 1
    resources:
        mem=low_memory_job,
    log:
        f"{logdir}GroupRefs_" "{sample}.log",
    params:
        script=srcdir("scripts/match_ref/group_refs.py") if config["use-conda"] is True and config["use-singularity"] is False else "/scripts/match_ref/group_refs.py",
    shell:
        """
        python {params.script} "{input.ref}" "{input.stats}" {output.groupedrefs} {output.groupedstats} {wildcards.sample} >> {log} 2>&1
        """


rule filter_gff:
    input:
        refdata=rules.group_and_rename_refs.output.groupedstats,
        gff=lambda wc: file if (file := SAMPLES[wc.sample]["FEATURES"]) else "",
    output:
        gff=f"{datadir}{matchref}" "{sample}_feats.gff",
        groupedstats=temp(f"{datadir}{matchref}" "{sample}_data1.csv"),
    threads: 1
    resources:
        mem=low_memory_job,
    log:
        f"{logdir}FilterGFF_" "{sample}.log",
    conda:
        f"{conda_envs}Scripts.yaml"
    container:
        f"{config['container_cache']}/viroconstrictor_scripts_{get_hash('Scripts')}.sif"
    params:
        script=srcdir("scripts/match_ref/filter_gff.py") if config["use-conda"] is True and config["use-singularity"] is False else "/scripts/match_ref/filter_gff.py",
    shell:
        """
        python {params.script} {input.refdata} {input.gff} {output.gff} {output.groupedstats} >> {log} 2>&1
        """

ruleorder: filter_gff > touch_gff

rule touch_gff:
    input:
        refdata=rules.group_and_rename_refs.output.groupedstats,
    output:
        gff=f"{datadir}{matchref}" "{sample}_feats.gff",
        groupedstats=temp(f"{datadir}{matchref}" "{sample}_data1.csv"),
    threads:1
    resources:
        mem=low_memory_job,
    shell:
        """
        cp {input.refdata} {output.groupedstats}
        touch {output.gff}
        """

rule filter_fasta2bed:
    input:
        ref=rules.group_and_rename_refs.output.groupedrefs,
        refdata=rules.filter_gff.output.groupedstats,
        prm=lambda wc: "" if SAMPLES[wc.sample]["PRIMERS"].endswith(".bed") else SAMPLES[wc.sample]["PRIMERS"],
    output:
        bed=f"{datadir}{matchref}" "{sample}_primers.bed",
        groupedstats=temp(f"{datadir}{matchref}" "{sample}_data2.csv"),
    conda:
        f"{conda_envs}Clean.yaml"
    container:
        f"{config['container_cache']}/viroconstrictor_clean_{get_hash('Clean')}.sif"
    threads: 1
    resources:
        mem_mb=low_memory_job,
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
        groupedstats=temp(f"{datadir}{matchref}" "{sample}_data2.csv"),
    threads: 1
    resources:
        mem_mb=low_memory_job,
    log:
        f"{logdir}FilterBed_" "{sample}.log",
    params:
        script=srcdir("scripts/match_ref/filter_bed.py") if config["use-conda"] is True and config["use-singularity"] is False else "/scripts/match_ref/filter_bed.py",
    conda:
        f"{conda_envs}Scripts.yaml"
    container:
        f"{config['container_cache']}/viroconstrictor_scripts_{get_hash('Scripts')}.sif"
    shell:
        """
        python {params.script} {input.prm} {input.refdata} {output.bed} {output.groupedstats}
        """


rule touch_primers:
    input:
        refdata=rules.filter_gff.output.groupedstats,
    output:
        bed=f"{datadir}{matchref}" "{sample}_primers.bed",
        groupedstats=temp(f"{datadir}{matchref}" "{sample}_data2.csv"),
    threads: 1
    resources:
        mem_mb=low_memory_job,
    shell:
        """
        cp {input.refdata} {output.groupedstats}
        touch {output.bed}
        """


rule concat_frames:
    input:
        set(
            expand(
                f"{datadir}{matchref}" "{sample}_data2.csv",
                sample=p_space.dataframe["sample"],
            )
        ),
    output:
        f"{datadir}" "match_ref_results.pkl",
    threads: 1
    resources:
        mem=low_memory_job,
    run:
        import pandas as pd

        df = pd.concat([pd.read_csv(f, keep_default_na=False) for f in input])
        df.to_pickle(output[0])


onsuccess:
    ViroConstrictor.logging.log.info(
        f"{'='*20} [green] Finished Match-reference process [/green] {'='*20}"
    )
    ViroConstrictor.logging.log.info(
        "[green]Finalizing the results and continuing with the main analysis workflow[/green]"
    )
    return True


onerror:
    ViroConstrictor.logging.log.error(
        "[bold red]An error occurred during the ViroConstrictor match-reference process.[/bold red]"
    )
    ViroConstrictor.logging.log.error(
        "[bold red]Shutting down... Please check all the inputs and logfiles for any abnormalities and try again.[/bold red]"
    )
    return False
