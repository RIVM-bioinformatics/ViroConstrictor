import pprint
import os
import yaml
import sys
import logging

import pandas as pd
import numpy as np

from directories import *
from presets import get_preset_parameter
from containers import get_hash
from Bio import SeqIO
import AminoExtract
from snakemake.utils import Paramspace, min_version

# setup and delete the log handlers which are imported from other modules.
# We only want the log handler from snakemake, not from the imported modules.
# import snakemake afterwards to re-set the log handlers without handlers from the imported modules.
logger = logging.getLogger()
logger.handlers.clear()

# you have to explicitly import the logger again, without this snakemake will not write to the log file.
from snakemake import logger

# Elevate the log level of all output generated by the snakemake.logging module to CRITICAL in order to suppress it when snakemake is calling itself in a downstream process.
if "--snakefile" in sys.argv:
    logging.getLogger("snakemake.logging").setLevel(logging.CRITICAL)
import ViroConstrictor

min_version("7.15")

yaml.warnings({"YAMLLoadWarning": False})
shell.executable("/bin/bash")

SAMPLES = {}
with open(config["sample_sheet"]) as sample_sheet_file:
    SAMPLES = yaml.safe_load(sample_sheet_file)


def Get_Ref_header(reffile):
    return [record.id for record in SeqIO.parse(reffile, "fasta")]


def Get_AA_feats(df):
    records = samples_df.to_dict(orient="records")

    for rec in records:
        if rec["FEATURES"] != "NONE":
            AA_dict = AminoExtract.get_feature_name_attribute(
                input_gff=str(rec["FEATURES"]),
                input_seq=str(rec["REFERENCE"]),
                feature_type="all",
            )
            if AA_dict:
                for k, v in AA_dict.items():
                    if k == rec["RefID"]:
                        rec["AA_FEAT_NAMES"] = tuple(v)
            else:
                rec["AA_FEAT_NAMES"] = np.nan
        else:
            rec["AA_FEAT_NAMES"] = np.nan
    return pd.DataFrame.from_records(records)


samples_df = (
    pd.DataFrame(SAMPLES)
    .transpose()
    .reset_index()
    .rename(columns=dict(index="sample", VIRUS="Virus"))
)
samples_df["RefID"] = samples_df["REFERENCE"].apply(Get_Ref_header)
samples_df = samples_df.explode("RefID")
samples_df = Get_AA_feats(samples_df)
p_space = Paramspace(
    samples_df[["Virus", "RefID", "sample"]], filename_params=["sample"]
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


localrules:
    all,
    prepare_refs,
    prepare_gffs,
    concat_sequences,
    concat_boc,
    concat_tsv_coverages,
    concat_amplicon_cov,
    make_pickle,
    move_fastq,
    create_empty_primers,


def list_aa_result_outputs():
    aa_feats = []
    for x in samples_df.to_dict(orient="records"):
        Virus = x["Virus"]
        RefID = x["RefID"]
        Feats = x["AA_FEAT_NAMES"]
        if not isinstance(Feats, float):
            aa_feats.extend(
                [f"{res}Virus~{Virus}/RefID~{RefID}/{amino}{aa}.faa" for aa in Feats]
            )
    return list(set(aa_feats))


def construct_all_rule(p_space):
    multiqc = f"{res}multiqc.html"
    folders = expand(
        f"{res}{wc_folder}",
        zip,
        RefID=p_space.RefID,
        Virus=p_space.Virus,
    )
    aa_feat_files = list_aa_result_outputs()

    base_results_files = expand(
        "{folder}{file}",
        folder=folders,
        file=[
            "consensus.fasta",
            "mutations.tsv",
            "Width_of_coverage.tsv",
            "Amplicon_coverage.csv",
        ],
    )

    return [multiqc] + base_results_files + aa_feat_files


wildcard_constraints:
    # regular expression to match only alphanumeric characters, underscores, dashes, and dots. exclude '/' and only match the first part of the string.
    RefID="[\w\-\.\d]+",
    Virus="[\w\-\.\d]+",
    # regular expression to match only alphanumeric characters, underscores, dashes. exclude '/' and only match the first part of the string.
    sample="[\w\-\.\d]+",


rule all:
    input:
        construct_all_rule(p_space),


rule prepare_refs:
    input:
        lambda wc: SAMPLES[wc.sample]["REFERENCE"],
    output:
        f"{datadir}{wc_folder}" "{sample}_reference.fasta",
    resources: 
        mem_mb=low_memory_job,
    threads: 1
    log:
        f"{logdir}prepare_refs_" "{Virus}.{RefID}.{sample}.log",
    benchmark:
        f"{logdir}{bench}prepare_refs_" "{Virus}.{RefID}.{sample}.txt"
    conda:
        f"{conda_envs}Scripts.yaml"
    container:
        f"{config['container_cache']}/viroconstrictor_scripts_{get_hash('Scripts')}.sif"
    params:
        script=srcdir("scripts/prepare_refs.py") if config["use-conda"] is True and config["use-singularity"] is False else "/scripts/prepare_refs.py",
    shell:
        """
        python {params.script} {input} {output} {wildcards.RefID} > {log}
        """

rule prepare_primers:
    input:
        prm=lambda wc: "" if SAMPLES[wc.sample]["PRIMERS"].endswith(".bed") else SAMPLES[wc.sample]["PRIMERS"],
        ref=rules.prepare_refs.output,
    output:
        bed=f"{datadir}{wc_folder}{prim}" "{sample}_primers.bed",
    resources:
        mem_mb=low_memory_job,
    log:
        f"{logdir}prepare_primers_" "{Virus}.{RefID}.{sample}.log",
    benchmark:
        f"{logdir}{bench}prepare_primers_" "{Virus}.{RefID}.{sample}.txt"
    params:
        pr_mm_rate=lambda wc: SAMPLES[wc.sample]["PRIMER-MISMATCH-RATE"],
    conda:
        f"{conda_envs}Clean.yaml"
    container:
        f"{config['container_cache']}/viroconstrictor_clean_{get_hash('Clean')}.sif"
    shell:
        """
        python -m AmpliGone.fasta2bed \
            --primers {input.prm} \
            --reference {input.ref} \
            --output {output.bed} \
            --primer-mismatch-rate {params.pr_mm_rate} \
            --verbose > {log}
        """

ruleorder: prepare_primers > filter_primer_bed > create_empty_primers

rule filter_primer_bed:
    input:
        prm=lambda wc: SAMPLES[wc.sample]["PRIMERS"],
    output:
        bed=f"{datadir}{wc_folder}{prim}" "{sample}_primers.bed",
    resources:
        mem_mb=low_memory_job,
    log:
        f"{logdir}prepare_primers_" "{Virus}.{RefID}.{sample}.log",
    benchmark:
        f"{logdir}{bench}prepare_primers_" "{Virus}.{RefID}.{sample}.txt"
    params:
        script=srcdir("scripts/filter_bed_input.py") if config["use-conda"] is True and config["use-singularity"] is False else "/scripts/filter_bed_input.py",
    conda:
        f"{conda_envs}Scripts.yaml"
    container:
        f"{config['container_cache']}/viroconstrictor_scripts_{get_hash('Scripts')}.sif"
    shell:
        """
        python {params.script} {input.prm} {output.bed} {wildcards.RefID}
        """

rule create_empty_primers:
    output:
        bed=touch(f"{datadir}{wc_folder}{prim}" "{sample}_primers.bed"),
    resources:
        mem_mb=low_memory_job,
    log:
        f"{logdir}prepare_primers_" "{Virus}.{RefID}.{sample}.log",
    shell:
        """
        echo "Created empty primer bed file for NONE primers" > {log}
        """


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
        f"{conda_envs}Scripts.yaml"
    container:
        f"{config['container_cache']}/viroconstrictor_scripts_{get_hash('Scripts')}.sif"
    resources:
        mem_mb=low_memory_job,
    params:
        script=srcdir("scripts/extract_gff.py") if config["use-conda"] is True and config["use-singularity"] is False else "/scripts/extract_gff.py",
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
        f"{conda_envs}ORF_analysis.yaml"
    container:
        f"{config['container_cache']}/viroconstrictor_orf_analysis_{get_hash('ORF_analysis')}.sif"
    resources:
        mem_mb=medium_memory_job,
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


if config["platform"] in ["nanopore", "iontorrent"]:
    p1_mapping_settings = (
        "-ax sr" if config["platform"] == "iontorrent" else "-ax map-ont"
    )

    rule qc_raw:
        input:
            lambda wc: SAMPLES[wc.sample]["INPUTFILE"],
        output:
            html=f"{datadir}{qc_pre}" "{sample}_fastqc.html",
            zip=f"{datadir}{qc_pre}" "{sample}_fastqc.zip",
        conda:
            f"{conda_envs}Scripts.yaml"
        container:
            f"{config['container_cache']}/viroconstrictor_scripts_{get_hash('Scripts')}.sif"
        log:
            f"{logdir}QC_raw_data_" "{sample}.log",
        benchmark:
            f"{logdir}{bench}QC_raw_data_" "{sample}.txt"
        threads: config["threads"]["QC"]
        resources:
            mem_mb=low_memory_job,
        params:
            output_dir=f"{datadir}{qc_pre}",
            script=srcdir("wrappers/fastqc_wrapper.sh") if config["use-conda"] is True and config["use-singularity"] is False else "/wrappers/fastqc_wrapper.sh",
        shell:
            """
            bash {params.script} {input} {params.output_dir} {output.html} {output.zip} {log}
            """

    rule remove_adapters_p1:
        input:
            ref=rules.prepare_refs.output,
            fq=lambda wc: SAMPLES[wc.sample]["INPUTFILE"],
        output:
            bam=f"{datadir}{wc_folder}{cln}{raln}" "{sample}.bam",
            index=f"{datadir}{wc_folder}{cln}{raln}" "{sample}.bam.bai",
        conda:
            f"{conda_envs}Alignment.yaml"
        container:
            f"{config['container_cache']}/viroconstrictor_alignment_{get_hash('Alignment')}.sif"
        log:
            f"{logdir}RemoveAdapters_p1_" "{Virus}.{RefID}.{sample}.log",
        benchmark:
            f"{logdir}{bench}RemoveAdapters_p1_" "{Virus}.{RefID}.{sample}.txt"
        threads: config["threads"]["Alignments"]
        resources:
            mem_mb=medium_memory_job,
        params:
            mapthreads=config["threads"]["Alignments"] - 1,
            mapping_base_settings=p1_mapping_settings,
            mapping_additionalsettings=lambda wc: get_preset_parameter(
                preset_name=SAMPLES[wc.sample]["PRESET"],
                parameter_name="RawAlign_AdditionalSettings",
            ),
            BAMfilters=lambda wc: get_preset_parameter(
                preset_name=SAMPLES[wc.sample]["PRESET"],
                parameter_name="BaseBAMFilters",
            ),
        shell:
            """
            minimap2 {params.mapping_base_settings} {params.mapping_additionalsettings} -t {params.mapthreads} {input.ref} {input.fq} 2>> {log} |\
            samtools view -@ {threads} {params.BAMfilters} -uS 2>> {log} |\
            samtools sort -o {output.bam} >> {log} 2>&1
            samtools index {output.bam} >> {log} 2>&1
            """


if config["platform"] == "illumina":

    rule qc_raw:
        input:
            lambda wildcards: SAMPLES[wildcards.sample][wildcards.read],
        output:
            html=f"{datadir}{qc_pre}" "{sample}_{read}_fastqc.html",
            zip=f"{datadir}{qc_pre}" "{sample}_{read}_fastqc.zip",
        conda:
            f"{conda_envs}Scripts.yaml"
        container:
            f"{config['container_cache']}/viroconstrictor_scripts_{get_hash('Scripts')}.sif"
        log:
            f"{logdir}" "QC_raw_data_{sample}_{read}.log",
        benchmark:
            f"{logdir}{bench}" "QC_raw_data_{sample}_{read}.txt"
        threads: config["threads"]["QC"]
        resources:
            mem_mb=low_memory_job,
        params:
            output_dir=f"{datadir}{qc_pre}",
            script=srcdir("wrappers/fastqc_wrapper.sh") if config["use-conda"] is True and config["use-singularity"] is False else "/wrappers/fastqc_wrapper.sh",
        shell:
            """
            bash {params.script} {input} {params.output_dir} {output.html} {output.zip} {log}
            """

    rule remove_adapters_p1:
        input:
            ref=rules.prepare_refs.output,
            fq1=lambda wildcards: SAMPLES[wildcards.sample]["R1"],
            fq2=lambda wildcards: SAMPLES[wildcards.sample]["R2"],
        output:
            bam=f"{datadir}{wc_folder}{cln}{raln}" "{sample}.bam",
            index=f"{datadir}{wc_folder}{cln}{raln}" "{sample}.bam.bai",
        conda:
            f"{conda_envs}Alignment.yaml"
        container:
            f"{config['container_cache']}/viroconstrictor_alignment_{get_hash('Alignment')}.sif"
        log:
            f"{logdir}" "RemoveAdapters_p1_{Virus}.{RefID}.{sample}.log",
        benchmark:
            f"{logdir}{bench}" "RemoveAdapters_p1_{Virus}.{RefID}.{sample}.log"
        threads: config["threads"]["Alignments"]
        resources:
            mem_mb=medium_memory_job,
        params:
            mapthreads=config["threads"]["Alignments"] - 1,
            mapping_additionalsettings=lambda wc: get_preset_parameter(
                preset_name=SAMPLES[wc.sample]["PRESET"],
                parameter_name="RawAlign_AdditionalSettings",
            ),
            BAMfilters=lambda wc: get_preset_parameter(
                preset_name=SAMPLES[wc.sample]["PRESET"],
                parameter_name="BaseBAMFilters",
            ),
        shell:
            """
            minimap2 -ax sr {params.mapping_additionalsettings} -t {params.mapthreads} {input.ref} {input.fq1} {input.fq2} 2>> {log} |\
            samtools view -@ {threads} {params.BAMfilters} -uS 2>> {log} |\
            samtools sort -o {output.bam} >> {log} 2>&1
            samtools index {output.bam} >> {log} 2>&1
            """


rule remove_adapters_p2:
    input:
        rules.remove_adapters_p1.output.bam,
    output:
        f"{datadir}{wc_folder}{cln}{noad}" "{sample}.fastq",
    conda:
        f"{conda_envs}Scripts.yaml"
    container:
        f"{config['container_cache']}/viroconstrictor_scripts_{get_hash('Scripts')}.sif"
    threads: config["threads"]["AdapterRemoval"]
    resources:
        mem_mb=low_memory_job,
    params:
        script=srcdir("scripts/clipper.py") if config["use-conda"] is True and config["use-singularity"] is False else "/scripts/clipper.py",
        clipper_settings=lambda wc: get_preset_parameter(
            preset_name=SAMPLES[wc.sample]["PRESET"],
            parameter_name=f"ClipperSettings_{config['platform']}",
        ),
    shell:
        """
        python {params.script} --input {input} --output {output} {params.clipper_settings} --threads {threads}
        """


rule qc_filter:
    input:
        rules.remove_adapters_p2.output,
    output:
        fq=f"{datadir}{wc_folder}{cln}{qcfilt}" "{sample}.fastq",
        html=f"{datadir}{wc_folder}{cln}{qcfilt}{html}" "{sample}_fastp.html",
        json=f"{datadir}{wc_folder}{cln}{qcfilt}{json}" "{sample}_fastp.json",
    conda:
        f"{conda_envs}Clean.yaml"
    container:
        f"{config['container_cache']}/viroconstrictor_clean_{get_hash('Clean')}.sif"
    log:
        f"{logdir}QC_filter_" "{Virus}.{RefID}.{sample}.log",
    benchmark:
        f"{logdir}{bench}QC_filter_" "{Virus}.{RefID}.{sample}.txt"
    threads: config["threads"]["QC"]
    resources:
        mem_mb=low_memory_job,
    params:
        score=lambda wc: get_preset_parameter(
            preset_name=SAMPLES[wc.sample]["PRESET"],
            parameter_name=f"Fastp_PhredScore_cutoff_{config['platform']}",
        ),
        size=lambda wc: get_preset_parameter(
            preset_name=SAMPLES[wc.sample]["PRESET"],
            parameter_name=f"Fastp_WindowSize_{config['platform']}",
        ),
        length=lambda wc: get_preset_parameter(
            preset_name=SAMPLES[wc.sample]["PRESET"],
            parameter_name="Fastp_MinReadLength",
        ),
    shell:
        """
        fastp --thread {threads} -i {input} \
            -A -Q --cut_right \
            --cut_right_mean_quality {params.score} \
            --cut_right_window_size {params.size} \
            -l {params.length} -o {output.fq} \
            -h {output.html} -j {output.json} > {log} 2>&1
        """

def get_primers_output(wildcards):
    sample_primers = SAMPLES[wildcards.sample]["PRIMERS"]
    if sample_primers == "NONE":
        return ""
    elif sample_primers.endswith(".bed"):
        return rules.filter_primer_bed.output.bed
    else:
        return rules.prepare_primers_fasta.output.bed

rule ampligone:
    input:
        fq=rules.qc_filter.output.fq,
        pr=lambda wc: get_primers_output(wc),
        ref=rules.prepare_refs.output,
    output:
        fq=f"{datadir}{wc_folder}{cln}{prdir}" "{sample}.fastq",
        ep=f"{datadir}{wc_folder}{prim}" "{sample}_removedprimers.bed",
    conda:
        f"{conda_envs}Clean.yaml"
    container:
        f"{config['container_cache']}/viroconstrictor_clean_{get_hash('Clean')}.sif"
    log:
        f"{logdir}" "AmpliGone_{Virus}.{RefID}.{sample}.log",
    benchmark:
        f"{logdir}{bench}" "AmpliGone_{Virus}.{RefID}.{sample}.txt"
    threads: config["threads"]["PrimerRemoval"]
    resources:
        mem_mb=high_memory_job,
    params:
        amplicontype=config["amplicon_type"],
        primer_mismatch_rate=lambda wc: SAMPLES[wc.sample]["PRIMER-MISMATCH-RATE"],
        alignmentpreset=lambda wc: get_preset_parameter(
            preset_name=SAMPLES[wc.sample]["PRESET"],
            parameter_name=f"AmpliGone_AlignmentPreset_{config['platform']}",
        ),
        alignmentmatrix=lambda wc: get_preset_parameter(
            preset_name=SAMPLES[wc.sample]["PRESET"],
            parameter_name=f"AmpliGone_AlignmentMatrix_{config['platform']}",
        ),
        extrasettings=lambda wc: get_preset_parameter(
            preset_name=SAMPLES[wc.sample]["PRESET"],
            parameter_name=f"AmpliGone_ExtraSettings",
        )
    shell:
        """
        echo "Running AmpliGone with primer file: {input.pr}" > {log}
        AmpliGone \
            -i {input.fq} \
            -ref {input.ref} -pr {input.pr} \
            -o {output.fq} \
            -at {params.amplicontype} \
            --error-rate {params.primer_mismatch_rate} \
            --export-primers {output.ep} \
            {params.alignmentpreset} {params.alignmentmatrix} {params.extrasettings} \
            -to \
            -t {threads} >> {log} 2>&1
        """


ruleorder: ampligone > move_fastq


# If no primers are given (e.g. with illumina runs), this rule makes sure the fastq's end up in the right place
rule move_fastq:
    input:
        rules.qc_filter.output.fq,
    output:
        fq=f"{datadir}{wc_folder}{cln}{prdir}" "{sample}.fastq",
        ep=touch(f"{datadir}{wc_folder}{prim}" "{sample}_removedprimers.bed"),
        # pr=touch(f"{datadir}{wc_folder}{prim}" "{sample}_primers.bed"),
    resources:
        mem_mb=low_memory_job,
    shell:
        """
        cp {input} {output.fq}
        """


rule qc_clean:
    input:
        rules.qc_filter.output.fq,
    output:
        html=f"{datadir}{wc_folder}{qc_post}" "{sample}_fastqc.html",
        zip=f"{datadir}{wc_folder}{qc_post}" "{sample}_fastqc.zip",
    conda:
        f"{conda_envs}Scripts.yaml"
    container:
        f"{config['container_cache']}/viroconstrictor_scripts_{get_hash('Scripts')}.sif"
    log:
        f"{logdir}QC_clean_" "{Virus}.{RefID}.{sample}.log",
    benchmark:
        f"{logdir}{bench}QC_clean_" "{Virus}.{RefID}.{sample}.txt"
    threads: config["threads"]["QC"]
    resources:
        mem_mb=low_memory_job,
    params:
        outdir=f"{datadir}{wc_folder}{qc_post}",
        script=srcdir("wrappers/fastqc_wrapper.sh") if config["use-conda"] is True and config["use-singularity"] is False else "/wrappers/fastqc_wrapper.sh",
    shell:
        """
        bash {params.script} {input} {params.outdir} {output.html} {output.zip} {log}
        """


### Align the cleaned reads to the reference
def get_alignment_flags(_wildcards):
    # These should correspond to the alignment settings in AmpliGone
    if config["platform"] in ["illumina", "iontorrent"]:
        return "-ax sr"
    elif config["platform"] == "nanopore":
        return "-ax map-ont"


rule align_before_trueconsense:
    input:
        fq=f"{datadir}{wc_folder}{cln}{prdir}" "{sample}.fastq",  #rules.ampligone.output.fq,
        ref=rules.prepare_refs.output,
    output:
        bam=f"{datadir}{wc_folder}{aln}{bf}" "{sample}.bam",
        index=f"{datadir}{wc_folder}{aln}{bf}" "{sample}.bam.bai",
    conda:
        f"{conda_envs}Alignment.yaml"
    container:
        f"{config['container_cache']}/viroconstrictor_alignment_{get_hash('Alignment')}.sif"
    log:
        f"{logdir}Alignment_" "{Virus}.{RefID}.{sample}.log",
    benchmark:
        f"{logdir}{bench}Alignment_" "{Virus}.{RefID}.{sample}.txt"
    threads: config["threads"]["Alignments"]
    resources:
        mem_mb=medium_memory_job,
    params:
        mapthreads=config["threads"]["Alignments"] - 1,
        alignment_base_settings=get_alignment_flags,
        alignment_additional_settings=lambda wc: get_preset_parameter(
            preset_name=SAMPLES[wc.sample]["PRESET"],
            parameter_name=f"CleanAlign_AdditionalSettings_{config['platform']}",
        ),
        BAMfilters=lambda wc: get_preset_parameter(
            preset_name=SAMPLES[wc.sample]["PRESET"],
            parameter_name=f"BaseBAMFilters",
        ),
    shell:
        """
        minimap2 {params.alignment_base_settings} {params.alignment_additional_settings} -t {params.mapthreads} {input.ref} {input.fq} 2>> {log} |\
        samtools view -@ {threads} {params.BAMfilters} -uS 2>> {log} |\
        samtools sort -o {output.bam} >> {log} 2>&1
        samtools index {output.bam} >> {log} 2>&1
        """


rule trueconsense:
    input:
        bam=rules.align_before_trueconsense.output.bam,
        gff=f"{datadir}{wc_folder}{features}" "{sample}_features.gff",
        ref=rules.prepare_refs.output,
    output:
        cons=f"{datadir}{wc_folder}{cons}{seqs}" "{sample}.fa",
        cov=f"{datadir}{wc_folder}{cons}{covs}" "{sample}_coverage.tsv",
        vcf=f"{datadir}{wc_folder}{aln}{vf}" "{sample}.vcf",
        gff=f"{datadir}{wc_folder}{cons}{features}" "{sample}.gff",
    params:
        mincov=lambda wc: SAMPLES[wc.sample]["MIN-COVERAGE"],
    conda:
        f"{conda_envs}Consensus.yaml"
    container:
        f"{config['container_cache']}/viroconstrictor_consensus_{get_hash('Consensus')}.sif"
    log:
        f"{logdir}Consensus_" "{Virus}.{RefID}.{sample}.log",
    benchmark:
        f"{logdir}{bench}Trueconsense_" "{Virus}.{RefID}.{sample}.txt"
    threads: config["threads"]["Consensus"]
    resources:
        mem_mb=medium_memory_job,
    shell:
        """
        TrueConsense --input {input.bam} \
        --reference {input.ref} --features {input.gff} \
        --coverage-level {params.mincov} \
        --samplename {wildcards.sample} \
        --output {output.cons} \
        --variants {output.vcf} \
        --output-gff {output.gff} \
        --depth-of-coverage {output.cov} \
        --threads {threads} > {log} 2>&1
        """


def group_items(wildcards, folder, filename):
    filtered_virus = p_space.dataframe.loc[
        p_space.dataframe["Virus"] == wildcards.Virus
    ]
    filtered_refid = filtered_virus.loc[filtered_virus["RefID"] == wildcards.RefID]
    return [f"{folder}{item}{filename}" for item in list(filtered_refid["sample"])]


rule concat_sequences:
    input:  # we don't use group_items here as mincov is in the filename
        lambda wc: (
            f"{datadir}{wc_folder}{cons}{seqs}{sample}.fa"
            for sample in p_space.dataframe.loc[
                (p_space.Virus == wc.Virus) & (p_space.RefID == wc.RefID)
            ]["sample"]
        ),
    output:
        f"{res}{wc_folder}consensus.fasta",
    resources:
        mem_mb=low_memory_job,
    shell:
        "cat {input} >> {output}"


rule Translate_AminoAcids:
    input:
        seq=rules.trueconsense.output.cons,
        gff=rules.trueconsense.output.gff,
    output:
        f"{datadir}{wc_folder}{amino}" "{sample}/aa.faa",
    conda:
        f"{conda_envs}ORF_analysis.yaml"
    container:
        f"{config['container_cache']}/viroconstrictor_orf_analysis_{get_hash('ORF_analysis')}.sif"
    resources:
        mem_mb=low_memory_job,
    log:
        f"{logdir}Translate_AA_" "{Virus}.{RefID}.{sample}.log",
    benchmark:
        f"{logdir}{bench}Translate_AA_" "{Virus}.{RefID}.{sample}.txt"
    params:
        feature_type="all",
    shell:
        """
        AminoExtract -i {input.seq} -gff {input.gff} -o {output} -ft {params.feature_type} -n {wildcards.sample} --keep-gaps >> {log} 2>&1
        """


def group_aminoacids_inputs(wildcards):
    filtered_df = samples_df.loc[samples_df["AA_FEAT_NAMES"].notnull()]
    filtered_vir_list = list(filtered_df["Virus"].unique())

    struct = {}
    for i in filtered_vir_list:
        select_samples = list(
            samples_df.loc[samples_df["Virus"] == i]["sample"].unique()
        )
        # for x in select_samples:
        #     y = samples_df.loc[(samples_df["Virus"] == i) & (samples_df["sample"] == x)]["RefID"].unique()
        #     print(y)
        # select_refIDs = list(samples_df.loc[samples_df["Virus"] == i]["RefID"].unique())
        # print(select_refIDs)

        # create a dictionary of dictionaries for each virus, with 'i' as the primary key and sample as the secondary key having a list of refIDs as the value
        struct[i] = {
            sample: list(
                samples_df.loc[
                    (samples_df["Virus"] == i) & (samples_df["sample"] == sample)
                ]["RefID"].unique()
            )
            for sample in select_samples
        }
    file_list = []
    for virus, sample in struct.items():
        for sample, refid in sample.items():
            for ref in refid:
                file_list.append(
                    f"{datadir}Virus~{virus}/RefID~{ref}/{amino}{sample}/aa.faa"
                )
    return file_list


# this rule cannot and should not run in a separate environment/container as the sole purpose is to transfer data of the paramspace into something that can then be used in a later rule.
rule make_pickle:
    output:
        temp(f"{datadir}sampleinfo.pkl"),
    resources:
        mem_mb=low_memory_job,
    threads: 1
    params:
        space=samples_df[~samples_df["AA_FEAT_NAMES"].isnull()],
    run:
        import pandas as pd

        params.space.to_pickle(output[0], compression=None)


rule concat_aminoacids:
    input:
        files=lambda wildcards: group_aminoacids_inputs(wildcards),
        sampleinfo=rules.make_pickle.output,
    output:
        list_aa_result_outputs(),
    resources:
        mem_mb=low_memory_job,
    conda:
        f"{conda_envs}Scripts.yaml"
    container:
        f"{config['container_cache']}/viroconstrictor_scripts_{get_hash('Scripts')}.sif"
    threads: 1
    params:
        script=srcdir("scripts/group_aminoacids.py") if config["use-conda"] is True and config["use-singularity"] is False else "/scripts/group_aminoacids.py",
    shell:
        'python {params.script} "{input.files}" "{output}" {input.sampleinfo}'


rule vcf_to_tsv:
    input:
        vcf=rules.trueconsense.output.vcf,
    output:
        tsv=temp(f"{datadir}{wc_folder}{aln}{vf}" "{sample}.tsv"),
    conda:
        f"{conda_envs}Scripts.yaml"
    container:
        f"{config['container_cache']}/viroconstrictor_scripts_{get_hash('Scripts')}.sif"
    threads: config["threads"]["Index"]
    resources:
        mem_mb=low_memory_job,
    log:
        f"{logdir}" "vcf_to_tsv_{Virus}.{RefID}.{sample}.log",
    params:
        script=srcdir("scripts/vcf_to_tsv.py") if config["use-conda"] is True and config["use-singularity"] is False else "/scripts/vcf_to_tsv.py",
    shell:
        """
        python {params.script} {input.vcf} {output.tsv} {wildcards.sample} >> {log} 2>&1
        """


rule concat_tsv_coverages:
    input:
        lambda wildcards: group_items(
            wildcards, folder=f"{datadir}{wc_folder}{aln}{vf}", filename=".tsv"
        ),
    output:
        f"{res}{wc_folder}mutations.tsv",
    resources:
        mem_mb=low_memory_job,
    shell:
        """
        echo -e 'Sample\tReference_ID\tPosition\tReference_Base\tVariant_Base\tDepth' > {output}
        cat {input} >> {output}
        """


rule get_breadth_of_coverage:
    input:
        reference=rules.prepare_refs.output,
        coverage=rules.trueconsense.output.cov,
    output:
        temp(f"{datadir}{wc_folder}{boc}" "{sample}.tsv"),
    resources:
        mem_mb=low_memory_job,
    conda:
        f"{conda_envs}Scripts.yaml"
    container:
        f"{config['container_cache']}/viroconstrictor_scripts_{get_hash('Scripts')}.sif"
    params:
        script=srcdir("scripts/boc.py") if config["use-conda"] is True and config["use-singularity"] is False else "/scripts/boc.py",
    shell:
        """
        python {params.script} {input.reference} {wildcards.sample} {input.coverage} {output}
        """


rule concat_boc:
    input:
        lambda wildcards: group_items(
            wildcards, folder=f"{datadir}{wc_folder}{boc}", filename=".tsv"
        ),
    output:
        f"{res}{wc_folder}Width_of_coverage.tsv",
    resources:
        mem_mb=low_memory_job,
    shell:
        """
        echo -e "Sample_name\tWidth_at_mincov_1\tWidth_at_mincov_5\tWidth_at_mincov_10\tWidth_at_mincov_50\tWidth_at_mincov_100" > {output}
        cat {input} >> {output}
        """


rule calculate_amplicon_cov:
    input:
        pr=f"{datadir}{wc_folder}{prim}" "{sample}_primers.bed",
        cov=rules.trueconsense.output.cov,
    output:
        f"{datadir}{wc_folder}{prim}" "{sample}_ampliconcoverage.csv",
    log:
        f"{logdir}" "calculate_amplicon_cov_{Virus}.{RefID}.{sample}.log",
    benchmark:
        f"{logdir}{bench}" "calculate_amplicon_cov_{Virus}.{RefID}.{sample}.txt"
    resources:
        mem_mb=low_memory_job,
    conda:
        f"{conda_envs}Scripts.yaml"
    container:
        f"{config['container_cache']}/viroconstrictor_scripts_{get_hash('Scripts')}.sif"
    params:
        script=srcdir("scripts/amplicon_covs.py") if config["use-conda"] is True and config["use-singularity"] is False else "/scripts/amplicon_covs.py",
    shell:
        """
        python {params.script} \
        --primers {input.pr} \
        --coverages {input.cov} \
        --key {wildcards.sample} \
        --output {output} > {log} 2>&1
        """


rule concat_amplicon_cov:
    input:
        lambda wildcards: group_items(
            wildcards,
            folder=f"{datadir}{wc_folder}{prim}",
            filename="_ampliconcoverage.csv",
        ),
    output:
        f"{res}{wc_folder}Amplicon_coverage.csv",
    resources:
        mem_mb=low_memory_job,
    conda:
        f"{conda_envs}Scripts.yaml"
    container:
        f"{config['container_cache']}/viroconstrictor_scripts_{get_hash('Scripts')}.sif"	
    params:
        script=srcdir("scripts/concat_amplicon_covs.py") if config["use-conda"] is True and config["use-singularity"] is False else "/scripts/concat_amplicon_covs.py",
    shell:
        """
        python {params.script} --output {output} --input {input}
        """


def construct_MultiQC_input(_wildcards):
    if config["platform"] == "nanopore" or config["platform"] == "iontorrent":
        pre = expand(
            f"{datadir}{qc_pre}" "{sample}_fastqc.zip",
            sample=SAMPLES,
        )
    elif config["platform"] == "illumina":
        pre = expand(
            f"{datadir}{qc_pre}" "{sample}_{read}_fastqc.zip",
            sample=SAMPLES,
            read="R1 R2".split(),
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
        f"{conda_envs}Clean.yaml"
    container:
        f"{config['container_cache']}/viroconstrictor_clean_{get_hash('Clean')}.sif"
    log:
        f"{logdir}MultiQC_report.log",
    benchmark:
        f"{logdir}{bench}MultiQC_report.txt"
    resources:
        mem_mb=high_memory_job,
    params:
        conffile=srcdir("files/multiqc_config.yaml") if config["use-conda"] is True and config["use-singularity"] is False else "/files/multiqc_config.yaml",
        outdir=res,
    shell:
        """
        multiqc -d --force --config {params.conffile} -o {params.outdir} -n multiqc.html {input} > {log} 2>&1
        """


onsuccess:
    ViroConstrictor.logging.log.info(
        "[bold green]ViroConstrictor is finished with processing all the files in the given input directory.[/bold green]"
    )
    ViroConstrictor.logging.log.info(
        "[bold green]Generating reports and shutting down...[/bold green]"
    )
    return True


onerror:
    ViroConstrictor.logging.log.error(
        "[bold red]An error occurred and ViroConstrictor had to shut down.[/bold red]"
    )
    ViroConstrictor.logging.log.error(
        "[bold red]Please check the input and logfiles for any abnormalities and try again.[/bold red]"
    )
    return False
