import pprint
import os
import yaml
import sys
from directories import *
import shutil
from Bio import SeqIO
from snakemake.utils import Paramspace, min_version
import pandas as pd

min_version("6.0")

yaml.warnings({'YAMLLoadWarning': False})
shell.executable("/bin/bash")

SAMPLES = {}
with open(config["sample_sheet"]) as sample_sheet_file:
    SAMPLES = yaml.safe_load(sample_sheet_file)

# primers_extension = config["primer_file"].split('.')[-1]

# reffile = config["reference_file"]
# ref_basename = os.path.splitext(os.path.basename(reffile))[0]
#
# features_file = config["features_file"]

mincov = 30

def Get_Ref_header(reffile):
    return [record.id for record in SeqIO.parse(reffile, "fasta")]

def construct_paramspace(sampleinfo):
    space = []
    for key, val in sampleinfo.items():
    
        if val["MATCH-REF"] is False:
            for id in Get_Ref_header(val["REFERENCE"]):
                space.append({"Virus": val["VIRUS"], "RefID": id, "Sample": key})
    return Paramspace(pd.DataFrame.from_dict(space))

paramspace = construct_paramspace(SAMPLES)

## copy references to right locations
for sample, vals in SAMPLES.items():
    prims = vals['PRIMERS']
    feats = vals['FEATURES']
    ref = vals['REFERENCE']
    vir = vals['VIRUS']
    prims_ext = prims.split('.')[-1]

    for id in Get_Ref_header(ref):
        base = f"{datadir}/{vir}/{id}/{sample}.raw_align.bam"
        shutil.copy(ref, f"{base}reference.fasta")
        shutil.copy(feats, f"{base}features.gff")
        shutil.copy(prims, f"{base}primers.{prims_ext}")

def construct_all_rule(sampleinfo):
    files = set()
    files.add(f"{res}multiqc.html")

    for key, val in sampleinfo.items():
        if val["MATCH-REF"] is False:
            for id in Get_Ref_header(val["REFERENCE"]):
                files.add(f"{res}{val['VIRUS']}/{id}/consensus.fasta")
                files.add(f"{res}{val['VIRUS']}/{id}/mutations.tsv")
                files.add(f"{res}{val['VIRUS']}/{id}/Width_of_coverage.tsv")
                if val["PRIMERS"] is not None:
                    files.add(f"{res}{val['VIRUS']}/{id}/Amplicon_coverage.tsv")

    return list(files)

def construct_MultiQC_input(_wildcards):
    if config['platform'] == "nanopore" or config['platform'] == "iontorrent":
        pre = expand(
            f"{datadir}{qc_pre}{{param_struct}}_fastqc.zip", param_struct=paramspace.instance_patterns)
    elif config['platform'] == "illumina":
        pre = expand(
            f"{datadir}{qc_pre}{{param_struct}}_{read}_fastqc.zip",
            param_struct=paramspace.instance_patterns,
            read = "R1 R2".split()
            )
    else:
        raise ValueError(f"Platform {config['platform']} not recognised. Choose one of [illumina, nanopore, iontorrent].")
    post = expand(f"{datadir}{qc_post}{{param_struct}}_fastqc.zip", param_struct=paramspace.instance_patterns)

    return pre + post

def low_memory_job(wildcards, threads, attempt):
    if config['computing_execution'] == 'local':
        return min(attempt * threads * 1 * 1000, config['max_local_mem'])
    return attempt * threads * 1 * 1000

def medium_memory_job(wildcards, threads, attempt):
    if config['computing_execution'] == 'local':
        return min(attempt * threads * 2 * 1000, config['max_local_mem'])
    return attempt * threads * 2 * 1000

def high_memory_job(wildcards, threads, attempt):
    if config['computing_execution'] == 'local':
        return min(attempt * threads * 4 * 1000, config['max_local_mem'])
    return attempt * threads * 4 * 1000

localrules:
    all,
    prepare_refs

rule all:
    input: construct_all_rule(SAMPLES)

rule prepare_refs:
    input: lambda wildcards: SAMPLES[str(wildcards.sample).split('~')[1]]["REFERENCE"]
    output: 
        ref = f"{datadir}{refdir}{{target}}/{{refID}}/{{sample}}_reference.fasta",
        ref_indx = f"{datadir}{refdir}{{target}}/{{refID}}/{{sample}}_reference.fasta.fai",
    conda: f"{conda_envs}Alignment.yaml"
    threads: config['threads']['Index']
    resources:
        mem_mb = low_memory_job
    params:
        script=srcdir("scripts/extract_ref.py")
    shell:
        """
        python {params.script} {input} {output} {wildcards.refID}
        samtools faidx {output.ref} -o {output.ref_indx}
        """


#TODO: maybe redo this to make sure this also works when no primers are given
rule prepare_primers:
    input: 
        prm = lambda wildcards: SAMPLES[str(wildcards.sample).split('~')[1]]["PRIMERS"],
        ref = rules.prepare_refs.output.ref,
    output:
        bed = f"{datadir}{refdir}{{target}}/{{refID}}/{{sample}}_primers.bed",
    threads: 1
    resources:
        mem_mb= low_memory_job
    log: f"{logdir}prepare_primers_{{target}}.{{refID}}.{{sample}}.log"
    benchmark: f"{logdir}{bench}prepare_primers_{{target}}.{{refID}}.{{sample}}.txt"
    params:
        pr_mm_rate = lambda wildcards: SAMPLES[str(wildcards.sample).split('~')[1]]["PRIMER-MISMATCH-RATE"]
    conda: f"{conda_envs}Clean.yaml"
    shell:
        """
        if [[ {input.prm} == *.bed ]]:
            cp {input.prm} {output.bed}
        else:
            python -m AmpliGone.fasta2bed \
                --primers {input.prm} \
                --reference {input.ref} \
                --output {output.bed} \
                --primer-mismatch-rate {params.pr_mm_rate} > {log}
        """

rule prepare_gffs:
    input: 
        feats = lambda wildcards: SAMPLES[str(wildcards.sample).split('~')[1]]["FEATURES"],
        ref = rules.prepare_refs.output.ref,
    output:
        gff = f"{datadir}{refdir}{{target}}/{{refID}}/{{sample}}_features.gff",
        aa = f"{datadir}{refdir}{{target}}/{{refID}}/{{sample}}_features.aa.fasta",
        nt = f"{datadir}{refdir}{{target}}/{{refID}}/{{sample}}_features.nt.fasta",
    log: f"{logdir}prepare_gffs_{{target}}.{{refID}}.{{sample}}.log"
    benchmark: f"{logdir}{bench}prepare_gffs_{{target}}.{{refID}}.{{sample}}.txt"
    threads: config["threads"]["Index"]
    conda:
        f"{conda_envs}ORF_analysis.yaml"
    resources:
        mem_mb= medium_memory_job
    params:
        script=srcdir("scripts/extract_gff.py"),
        prodigal_method = "meta",
        prodigal_outformat= "gff",
    shell:
        """
        if [[ {input.feats} == "NONE" ]]:
            prodigal -q -i {input.ref} \
                -a {output.aa} \
                -d {output.nt} \
                -o {output.gff} \
                -p {params.prodigal_method} \
                -f {params.prodigal_outformat} > {log} 2>&1
        else:
            python {params.script} {input.feats} {output.gff} {wildcards.refID}
        """

if config["platform"] in ["nanopore", "iontorrent"]:
    p1_mapping_settings = "-ax sr" if config["platform"] == "iontorrent" else "-ax map-ont"

    rule qc_raw:
        input: lambda wildcards: SAMPLES[wildcards.sample]
        output:
            html = f"{datadir}{qc_pre}""{sample}_fastqc.html",
            zip  = f"{datadir}{qc_pre}""{sample}_fastqc.zip"
        conda: f"{conda_envs}Clean.yaml"
        log: f"{logdir}""QC_raw_data_{sample}.log"
        benchmark: f"{logdir}{bench}""QC_raw_data_{sample}.txt"
        threads: config['threads']['QC']
        resources: mem_mb = low_memory_job
        params:
            output_dir = f"{datadir}{qc_pre}",
            script = srcdir("scripts/fastqc_wrapper.sh")
        shell:
            """
            bash {params.script} {input} {params.output_dir} {output.html} {output.zip} {log}
            """

    rule remove_adapters_p1:
        input:
            ref = rules.prepare_ref.output.ref,
            fq  = lambda wildcards: SAMPLES[wildcards.sample]
        output:
            bam   = f"{datadir}{cln}{raln}""{sample}.bam",
            index = f"{datadir}{cln}{raln}""{sample}.bam.bai"
        conda: f"{conda_envs}Alignment.yaml"
        log: f"{logdir}""RemoveAdapters_p1_{sample}.log"
        benchmark: f"{logdir}{bench}"+ "RemoveAdapters_p1_{sample}.txt"
        threads: config['threads']['Alignments']
        resources: mem_mb = medium_memory_job
        params:
            mapthreads = config['threads']['Alignments'] - 1,
            filters = config["runparams"]["alignmentfilters"],
            mapping_settings = p1_mapping_settings
        shell:
            """
            minimap2 {params.mapping_settings} -t {params.mapthreads} {input.ref} {input.fq} 2>> {log} |\
            samtools view -@ {threads} {params.filters} -uS 2>> {log} |\
            samtools sort -o {output.bam} >> {log} 2>&1
            samtools index {output.bam} >> {log} 2>&1
            """

if config["platform"] == "illumina":
    rule qc_raw:
        input: lambda wildcards: SAMPLES[wildcards.sample][wildcards.read]
        output:
            html = f"{datadir}{qc_pre}""{sample}_{read}_fastqc.html",
            zip  = f"{datadir}{qc_pre}""{sample}_{read}_fastqc.zip",
        conda: f"{conda_envs}Clean.yaml"
        log: f"{logdir}""QC_raw_data_{sample}_{read}.log"
        benchmark: f"{logdir}{bench}""QC_raw_data_{sample}_{read}.txt"
        threads: config['threads']['QC']
        resources: mem_mb = low_memory_job
        params:
            output_dir = f"{datadir}{qc_pre}",
            script = srcdir("scripts/fastqc_wrapper.sh"),
        shell:
            """
            bash {params.script} {input} {params.output_dir} {output.html} {output.zip} {log}
            """

    rule remove_adapters_p1:
        input:
            ref = rules.prepare_ref.output.ref,
            fq  = lambda wildcards: (SAMPLES[wildcards.sample][i]
                                for i in ("R1", "R2")
                                )
        output:
            bam   = f"{datadir}{cln}{raln}""{sample}.bam",
            index = f"{datadir}{cln}{raln}""{sample}.bam.bai"
        conda: f"{conda_envs}Alignment.yaml"
        log: f"{logdir}""RemoveAdapters_p1_{sample}.log"
        benchmark: f"{logdir}{bench}""RemoveAdapters_p1_{sample}.txt"
        threads: config['threads']['Alignments']
        resources: mem_mb = medium_memory_job
        params:
            mapthreads = config['threads']['Alignments'] - 1,
            filters = config["runparams"]["alignmentfilters"]
        shell:
            """
            minimap2 -ax sr -t {params.mapthreads} {input.ref} {input.fq[0]:q} {input.fq[1]:q} 2>> {log} |\
            samtools view -@ {threads} {params.filters} -uS 2>> {log} |\
            samtools sort -o {output.bam} >> {log} 2>&1
            samtools index {output.bam} >> {log} 2>&1
            """

rule remove_adapters_p2:
    input: rules.remove_adapters_p1.output.bam
    output: f"{datadir}{cln}{noad}""{sample}.fastq"
    conda: f"{conda_envs}Clean.yaml"
    threads: config['threads']['AdapterRemoval']
    resources: mem_mb = low_memory_job
    params: script = srcdir('scripts/clipper.py')
    shell:
        """
        python {params.script} --input {input} --output {output} --threads {threads}
        """

rule qc_filter:
    input: rules.remove_adapters_p2.output
    output:
        fq = f"{datadir}{cln}{qcfilt}""{sample}.fastq",
        html = f"{datadir}{cln}{qcfilt}{html}""{sample}.fastp.html",
        json = f"{datadir}{cln}{qcfilt}{json}""{sample}.fastp.json"
    conda: f"{conda_envs}Clean.yaml"
    log: f"{logdir}""Cleanup_{sample}.log"
    benchmark: f"{logdir}{bench}""Cleanup_{sample}.txt"
    threads: config['threads']['QC']
    resources: mem_mb = low_memory_job
    params:
        score = config['runparams'][f"qc_filter_{config['platform']}"],
        size = config['runparams'][f"qc_window_{config['platform']}"],
        length = config['runparams']['qc_min_readlength']
    shell:
        """
        fastp --thread {threads} -i {input} \
            -A -Q --cut_right \
            --cut_right_mean_quality {params.score} \
            --cut_right_window_size {params.size} \
            -l {params.length} -o {output.fq} \
            -h {output.html} -j {output.json} > {log} 2>&1
        """

ruleorder: ampligone > move_fastq

rule ampligone:
    input:
        fq = rules.qc_filter.output.fq,
        pr = f"{datadir}{prim}primers.bed",
        ref = rules.prepare_ref.output.ref
    output:
        fq = f"{datadir}{cln}{prdir}""{sample}.fastq",
        ep = f"{datadir}{prim}""{sample}_removedprimers.bed"
    conda: f"{conda_envs}Clean.yaml"
    log: f"{logdir}""AmpliGone_{sample}.log"
    benchmark: f"{logdir}{bench}""AmpliGone_{sample}.txt"
    threads: config['threads']['PrimerRemoval']
    resources: mem_mb = high_memory_job
    params:
        amplicontype = config["amplicon_type"],
        pr_mm_rate = config["primer_mismatch_rate"]
    shell:
        """
        echo {input.pr} > {log}
        AmpliGone \
            -i {input.fq} \
            -ref {input.ref} -pr {input.pr} \
            -o {output.fq} \
            -at {params.amplicontype} \
            --error-rate {params.pr_mm_rate} \
            --export-primers {output.ep} \
            -to \
            -t {threads} >> {log} 2>&1
        """

rule move_fastq:
    input: rules.qc_filter.output.fq
    output: f"{datadir}{cln}{prdir}""{sample}.fastq"
    threads: 1
    resources: mem_mb = low_memory_job
    shell:
        """
        cp {input} {output}
        """

rule qc_clean:
    input: f"{datadir}{cln}{prdir}""{sample}.fastq"
    output:
        html = f"{datadir}{qc_post}""{sample}_fastqc.html",
        zip = f"{datadir}{qc_post}""{sample}_fastqc.zip"
    conda: f"{conda_envs}Clean.yaml"
    log: f"{logdir}""QC_clean_data_{sample}.log"
    benchmark: f"{logdir}{bench}""QC_clean_data_{sample}.txt"
    threads: config['threads']['QC']
    resources: mem_mb = low_memory_job
    params: outdir = f"{datadir}{qc_post}"
    shell:
        """
        if [ -s "{input}" ]; then
            fastqc -t {threads} --quiet --outdir {params.outdir} {input} > {log} 2>&1
        else
            touch {output.html}
            touch {output.zip}
        fi
        """


### Align the cleaned reads to the reference
def get_alignment_flags(_wildcards):
    # These should correspond to the alignment settings in AmpliGone
    if config["platform"] in ["illumina", "iontorrent"]:
        return "-ax sr"
    elif config["platform"] == "nanopore":
        return "-ax map-ont -E2,0 -O8,24 -A4 -B4"

rule align_before_trueconsense:
    input:
        fq = f"{datadir}{cln}{prdir}""{sample}.fastq",
        ref = rules.prepare_ref.output.ref
    output:
        bam = f"{datadir}{aln}{bf}""{sample}.bam",
        index = f"{datadir}{aln}{bf}""{sample}.bam.bai"
    conda: f"{conda_envs}Alignment.yaml"
    log: f"{logdir}""Alignment_{sample}.log"
    benchmark: f"{logdir}{bench}""Alignment_{sample}.txt"
    threads: config['threads']['Alignments']
    resources: mem_mb = medium_memory_job
    params:
        mapthreads = config['threads']['Alignments'] - 1,
        filters = config["runparams"]["alignmentfilters"],
        alignment_flags = get_alignment_flags,
    shell:
        """
        minimap2 {params.alignment_flags} -t {params.mapthreads} {input.ref} {input.fq} 2>> {log} |\
        samtools view -@ {threads} {params.filters} -uS 2>> {log} |\
        samtools sort -o {output.bam} >> {log} 2>&1
        samtools index {output.bam} >> {log} 2>&1
        """

rule trueconsense:
    input:
        bam = rules.align_before_trueconsense.output.bam,
        gff = rules.prepare_features_file.output.gff,
        ref = rules.prepare_ref.output.ref
    output:
        cons = f"{datadir}{cons}{seqs}""{sample}"f"_cov_ge_{mincov}.fa",
        cov = f"{datadir}{cons}{covs}""{sample}_coverage.tsv",
        vcf = f"{datadir}{aln}{vf}""{sample}"f"_cov_ge_{mincov}.vcf",
        gff = f"{datadir}{cons}{features}""{sample}"f"_cov_ge_{mincov}.gff",
    params:
        mincov = mincov,
        outdir = f"{datadir}{cons}{seqs}",
        vcfdir = f"{datadir}{aln}{vf}",
        gffdir = f"{datadir}{cons}{features}"
    conda: f"{conda_envs}Consensus.yaml"
    log: f"{logdir}""Consensus_{sample}.log"
    benchmark: f"{logdir}{bench}""Consensus_{sample}.txt"
    threads: config['threads']['Consensus']
    resources: mem_mb = medium_memory_job
    shell:
        """
        TrueConsense --input {input.bam} \
        --reference {input.ref} --features {input.gff} \
        --coverage-levels {params.mincov} \
        --samplename {wildcards.sample} \
        --output {params.outdir} \
        --variants {params.vcfdir} \
        --output-gff {params.gffdir} \
        --depth-of-coverage {output.cov} \
        --threads {threads}
        """

rule concat_sequences:
    input:
        expand(
            f"{datadir}{cons}{seqs}""{sample}_cov_ge_"f"{mincov}.fa",
            sample = SAMPLES,
            )
    output: f"{res}consensus.fasta",
    threads: 1
    resources: mem_mb = low_memory_job
    shell: "cat {input} >> {output}"


rule vcf_to_tsv:
    input: vcf = f"{datadir}{aln}{vf}""{sample}"f"_cov_ge_{mincov}.vcf",
    output: tsv = temp(f"{datadir}{aln}{vf}""{sample}.tsv"),
    conda: f"{conda_envs}Mutations.yaml"
    threads: config['threads']['Index']
    resources: mem_mb = low_memory_job
    log: f"{logdir}""vcf_to_tsv_{sample}.log"
    shell:
        """
        bcftools query {input.vcf} -f '{wildcards.sample}\t%CHROM\t%POS\t%REF\t%ALT\t%DP\n' -e 'ALT="N"' > {output.tsv} 2>> {log}
        """

rule concat_tsv_coverages:
    input: tsvs = expand(f"{datadir}{aln}{vf}""{sample}.tsv", sample = SAMPLES),
    output: tsv = f"{res}mutations.tsv",
    log: f"{logdir}concat_tsv.log"
    threads: 1
    resources: mem_mb = low_memory_job
    run:
        shell("echo -e 'Sample\tReference_Chromosome\tPosition\tReference\tAlternative\tDepth' > {output.tsv} 2> {log}")
        shell("cat {input.tsvs} >> {output.tsv} 2>> {log}")

rule get_breadth_of_coverage:
    input:
        reference = rules.prepare_ref.output.ref,
        coverage = rules.trueconsense.output.cov,
    output: temp(f"{datadir}{boc}""{sample}.tsv")
    conda: f"{conda_envs}Consensus.yaml"
    threads: 1
    resources: mem_mb = low_memory_job
    params: script = srcdir("scripts/boc.py")
    shell:
        """
        python {params.script} {input.reference} {wildcards.sample} {input.coverage} {output}
        """

rule concat_boc:
    input: pct = expand(f"{datadir}{boc}""{sample}.tsv", sample = SAMPLES)
    output: f"{res}Width_of_coverage.tsv"
    threads: 1
    resources: mem_mb = low_memory_job
    shell:
        """
        echo -e "Sample_name\tWidth_at_mincov_1\tWidth_at_mincov_5\tWidth_at_mincov_10\tWidth_at_mincov_50\tWidth_at_mincov_100" > {output}
        cat {input} >> {output}
        """

rule calculate_amplicon_cov:
    input:
        pr = rules.ampligone.output.ep,
        cov = rules.trueconsense.output.cov
    output: f"{datadir}{prim}""{sample}_ampliconcoverage.csv"
    threads: 1
    resources: mem_mb = low_memory_job
    conda: f"{conda_envs}Clean.yaml"
    params: script = srcdir("scripts/amplicon_covs.py")
    shell:
        """
        python {params.script} \
        --primers {input.pr} \
        --coverages {input.cov} \
        --key {wildcards.sample} \
        --output {output}
        """

rule concat_amplicon_cov:
    input: expand(f"{datadir}{prim}""{sample}_ampliconcoverage.csv", sample = SAMPLES)
    output: f"{res}Amplicon_coverage.csv"
    threads: 1
    resources: mem_mb = low_memory_job
    conda: f"{conda_envs}Clean.yaml"
    params: script = srcdir("scripts/concat_amplicon_covs.py")
    shell:
        """
        python {params.script} --output {output} --input {input}
        """

rule multiqc_report:
    input: construct_MultiQC_input
    output:
        f"{res}multiqc.html",
        expand(f"{res}{mqc_data}""multiqc_{program}.txt", program = "fastqc")
    conda: f"{conda_envs}Clean.yaml"
    log: f"{logdir}MultiQC_report.log"
    benchmark: f"{logdir}{bench}MultiQC_report.txt"
    threads: 1
    resources: mem_mb = medium_memory_job
    params:
        conffile = srcdir('files/multiqc_config.yaml'),
        outdir = res
    shell:
        """
        multiqc -d --force --config {params.conffile} -o {params.outdir} -n multiqc.html {input} > {log} 2>&1
        """

onsuccess:
    print("""
    ViroConstrictor is finished with processing all the files in the given input directory.

    Generating reports and shutting down...
    """)
    return True

onerror:
    print("""
    An error occurred and ViroConstrictor had to shut down.
    Please check the input and logfiles for any abnormalities and try again.
    """)
    return False
