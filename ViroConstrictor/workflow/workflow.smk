import pprint
import os
import yaml
import sys
from directories import *
import snakemake

snakemake.utils.min_version("6.0")

yaml.warnings({'YAMLLoadWarning': False})
shell.executable("/bin/bash")

SAMPLES = {}
with open(config["sample_sheet"]) as sample_sheet_file:
    SAMPLES = yaml.load(sample_sheet_file)

primerfile = config["primer_file"]

if primerfile == "NONE":
    primerfile = srcdir("files/empty.primers")

reffile = config["reference_file"]
ref_basename = os.path.splitext(os.path.basename(reffile))[0]

features_file = config["features_file"]

mincoverages = [1,5,10,50,100]

if config["primer_file"] == "NONE":
    rule all:
        input:
            f"{res}multiqc.html",
            expand("{p}concat_cov_ge_{cov}.fasta",
                p = res + seqs,
                cov = mincoverages),
            expand("{p}concat_mutations_cov_ge_{cov}.tsv",
                p = res + muts,
                cov = mincoverages),
            f"{res}Width_of_coverage.tsv",
else:
    rule all:
        input:
            f"{res}multiqc.html",
            expand("{p}concat_cov_ge_{cov}.fasta",
                p = res + seqs,
                cov = mincoverages),
            expand("{p}concat_mutations_cov_ge_{cov}.tsv",
                p = res + muts,
                cov = mincoverages),
            f"{res}Width_of_coverage.tsv",
            f"{res}Amplicon_coverage.csv"


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

#TODO: Check if this can be done more elegantly than by wrapping everything in an if/else statement
if config["primer_file"] != "NONE":
    rule Prepare_ref_and_primers:
        input:
            ref = reffile,
            prm = primerfile
        output: 
            ref = f"{datadir + refdir + ref_basename}.fasta",
            prm = f"{datadir + prim}" + "primers.fasta",
            refindex = f"{datadir + refdir + ref_basename}.fasta.fai"
        conda:
            f"{conda_envs}Alignment.yaml"
        threads: config['threads']['Index']
        resources:
            mem_mb = low_memory_job
        shell:
            """
            cat {input.ref} | seqkit replace -p "\-" -s -r "N" > {output.ref}
            cat {input.prm} | seqkit replace -p "\-" -s -r "N" > {output.prm}
            samtools faidx {output.ref} -o {output.refindex}
            """
else:
    rule Prepare_ref_and_primers:
        input:
            ref = reffile,
            prm = primerfile
        output: 
            ref = f"{datadir + refdir + ref_basename}.fasta",
            prm = temp(f"{datadir + prim}" + "primers.fasta"),
            refindex = f"{datadir + refdir + ref_basename}.fasta.fai"
        conda:
            f"{conda_envs}Alignment.yaml"
        threads: config['threads']['Index']
        resources:
            mem_mb = low_memory_job
        shell:
            """
            cat {input.ref} | seqkit replace -p "\-" -s -r "N" > {output.ref}
            cat {input.prm} | seqkit replace -p "\-" -s -r "N" > {output.prm}
            samtools faidx {output.ref} -o {output.refindex}
            """

if config["features_file"] != "NONE":
    rule Prepare_features_file:
        input: features_file
        output: 
            gff = temp(f"{datadir + features}reference_features.gff")
        threads: config['threads']['Index']
        resources:
            mem_mb = low_memory_job
        shell:
            """
            cat {input} >> {output}
            """
else:
    rule Prepare_features_file:
        input: rules.Prepare_ref_and_primers.output.ref
        output:
            AA  =   temp(f"{datadir + refdir + ref_basename}_ORF_AA.fa"),
            NT  =   temp(f"{datadir + refdir + ref_basename}_ORF_NT.fa"),
            gff =   temp(f"{datadir + refdir + ref_basename}_annotation.gff")
        conda:
            f"{conda_envs}ORF_analysis.yaml"
        threads: config['threads']['Index']
        resources:
            mem_mb = medium_memory_job
        log:
            f"{logdir}ORF_Analysis.log"
        params:
            procedure = "meta",
            output_format = "gff"
        shell:
            """
            prodigal -q -i {input} \
            -a {output.AA} \
            -d {output.NT} \
            -o {output.gff} \
            -p {params.procedure} \
            -f {params.output_format} > {log} 2>&1
            """



if config["platform"] == "illumina":
    rule QC_raw:
        input: lambda wildcards: SAMPLES[wildcards.sample][wildcards.read]
        output:
            html    =   f"{datadir + qc_pre}" + "{sample}_{read}_fastqc.html",
            zip     =   f"{datadir + qc_pre}" + "{sample}_{read}_fastqc.zip"
        conda:
            f"{conda_envs}Clean.yaml"
        log:
            f"{logdir}" + "QC_raw_data_{sample}_{read}.log"
        benchmark:
            f"{logdir + bench}" + "QC_raw_data_{sample}_{read}.txt"
        threads: config['threads']['QC']
        resources:
            mem_mb = low_memory_job
        params:
            output_dir  =   f"{datadir + qc_pre}",
            script = srcdir("scripts/fastqc_wrapper.sh")
        shell:
            """
            bash {params.script} {input} {params.output_dir} {output.html} {output.zip} {log}
            """

    rule RemoveAdapters_p1:
        input: 
            ref = rules.Prepare_ref_and_primers.output.ref,
            fq  = lambda wildcards: (SAMPLES[wildcards.sample][i]
                                for i in ("R1", "R2")
                                )
        output:
            bam     = f"{datadir + cln + raln}" + "{sample}.bam",
            index   = f"{datadir + cln + raln}" + "{sample}.bam.bai"
        conda:
            f"{conda_envs}Alignment.yaml"
        log:
            f"{logdir}" + "RemoveAdapters_p1_{sample}.log"
        benchmark:
            f"{logdir + bench}"+ "RemoveAdapters_p1_{sample}.txt"
        threads: config['threads']['Alignments']
        resources:
            mem_mb = medium_memory_job
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
    
    rule RemoveAdapters_p2:
        input: rules.RemoveAdapters_p1.output.bam
        output: f"{datadir + cln + noad}" + "{sample}.fastq"
        conda:
            f"{conda_envs}Clean.yaml"
        threads: config['threads']['AdapterRemoval']
        resources:
            mem_mb = low_memory_job
        params:
            script = srcdir('scripts/clipper.py')
        shell: 
            """
            python {params.script} --input {input} --output {output} --threads {threads}
            """

    rule QC_filter:
        input: rules.RemoveAdapters_p2.output
        output: 
            fq = f"{datadir + cln + qcfilt}" + "{sample}.fastq",
            html = f"{datadir + cln + qcfilt + html}" + "{sample}.fastp.html",
            json = f"{datadir + cln + qcfilt + json}" + "{sample}.fastp.json"
        conda: 
            f"{conda_envs}Clean.yaml"
        log:
            f"{logdir}" + "Cleanup_{sample}.log"
        benchmark:
            f"{logdir + bench}" + "Cleanup_{sample}.txt"
        threads: config['threads']['QC']
        resources:
            mem_mb = low_memory_job
        params:
            score = config['runparams']['qc_filter_illumina'],
            size = config['runparams']['qc_window_illumina'],
            length = config['runparams']['qc_min_readlength']
        shell:
            """
            fastp --thread {threads} -i {input} \
            --cut_right -M {params.score} -W {params.size} -l {params.length} \
            -o {output.fq} -h {output.html} -j {output.json} > {log} 2>&1
            """

if config["platform"] == "nanopore":
    rule QC_raw:
        input: lambda wildcards: SAMPLES[wildcards.sample]
        output:
            html    =   f"{datadir + qc_pre}" + "{sample}_fastqc.html",
            zip     =   f"{datadir + qc_pre}" + "{sample}_fastqc.zip"
        conda:
            f"{conda_envs}Clean.yaml"
        log:
            f"{logdir}" + "QC_raw_data_{sample}.log"
        benchmark:
            f"{logdir + bench}" + "QC_raw_data_{sample}.txt"
        threads: config['threads']['QC']
        resources:
            mem_mb = low_memory_job
        params:
            output_dir  =   f"{datadir + qc_pre}",
            script = srcdir("scripts/fastqc_wrapper.sh")
        shell:
            """
            bash {params.script} {input} {params.output_dir} {output.html} {output.zip} {log}
            """

    rule RemoveAdapters_p1:
        input: 
            ref = rules.Prepare_ref_and_primers.output.ref,
            fq  = lambda wildcards: SAMPLES[wildcards.sample]
        output:
            bam     = f"{datadir + cln + raln}" + "{sample}.bam",
            index   = f"{datadir + cln + raln}" + "{sample}.bam.bai"
        conda:
            f"{conda_envs}Alignment.yaml"
        log:
            f"{logdir}" + "RemoveAdapters_p1_{sample}.log"
        benchmark:
            f"{logdir + bench}"+ "RemoveAdapters_p1_{sample}.txt"
        threads: config['threads']['Alignments']
        resources:
            mem_mb = medium_memory_job
        params:
            mapthreads = config['threads']['Alignments'] - 1,
            filters = config["runparams"]["alignmentfilters"]
        shell: 
            """
            minimap2 -ax map-ont -t {params.mapthreads} {input.ref} {input.fq} 2>> {log} |\
            samtools view -@ {threads} {params.filters} -uS 2>> {log} |\
            samtools sort -o {output.bam} >> {log} 2>&1
            samtools index {output.bam} >> {log} 2>&1
            """
    
    rule RemoveAdapters_p2:
        input: rules.RemoveAdapters_p1.output.bam
        output: f"{datadir + cln + noad}" + "{sample}.fastq"
        conda:
            f"{conda_envs}Clean.yaml"
        threads: config['threads']['AdapterRemoval']
        resources:
            mem_mb = low_memory_job
        params:
            script = srcdir('scripts/clipper.py')
        shell: 
            """
            python {params.script} --input {input} --output {output} --threads {threads}
            """

    rule QC_filter:
        input: rules.RemoveAdapters_p2.output
        output: 
            fq = f"{datadir + cln + qcfilt}" + "{sample}.fastq",
            html = f"{datadir + cln + qcfilt + html}" + "{sample}.fastp.html",
            json = f"{datadir + cln + qcfilt + json}" + "{sample}.fastp.json"
        conda: 
            f"{conda_envs}Clean.yaml"
        log:
            f"{logdir}" + "Cleanup_{sample}.log"
        benchmark:
            f"{logdir + bench}" + "Cleanup_{sample}.txt"
        threads: config['threads']['QC']
        resources:
            mem_mb = low_memory_job
        params:
            score = config['runparams']['qc_filter_nanopore'],
            size = config['runparams']['qc_window_nanopore'],
            length = config['runparams']['qc_min_readlength']
        shell:
            """
            fastp --thread {threads} -i {input} \
            --cut_right -M {params.score} -W {params.size} -l {params.length} \
            -o {output.fq} -h {output.html} -j {output.json} > {log} 2>&1
            """

if config["platform"] == "iontorrent":
    rule QC_raw:
        input: lambda wildcards: SAMPLES[wildcards.sample]
        output:
            html    =   f"{datadir + qc_pre}" + "{sample}_fastqc.html",
            zip     =   f"{datadir + qc_pre}" + "{sample}_fastqc.zip"
        conda:
            f"{conda_envs}Clean.yaml"
        log:
            f"{logdir}" + "QC_raw_data_{sample}.log"
        benchmark:
            f"{logdir + bench}" + "QC_raw_data_{sample}.txt"
        threads: config['threads']['QC']
        resources:
            mem_mb = low_memory_job
        params:
            output_dir  =   f"{datadir + qc_pre}",
            script = srcdir("scripts/fastqc_wrapper.sh")
        shell:
            """
            bash {params.script} {input} {params.output_dir} {output.html} {output.zip} {log}
            """

    rule RemoveAdapters_p1:
        input: 
            ref = rules.Prepare_ref_and_primers.output.ref,
            fq  = lambda wildcards: SAMPLES[wildcards.sample]
        output:
            bam     = f"{datadir + cln + raln}" + "{sample}.bam",
            index   = f"{datadir + cln + raln}" + "{sample}.bam.bai"
        conda:
            f"{conda_envs}Alignment.yaml"
        log:
            f"{logdir}" + "RemoveAdapters_p1_{sample}.log"
        benchmark:
            f"{logdir + bench}"+ "RemoveAdapters_p1_{sample}.txt"
        threads: config['threads']['Alignments']
        resources:
            mem_mb = medium_memory_job
        params:
            mapthreads = config['threads']['Alignments'] - 1,
            filters = config["runparams"]["alignmentfilters"]
        shell: 
            """
            minimap2 -ax sr -t {params.mapthreads} {input.ref} {input.fq} 2>> {log} |\
            samtools view -@ {threads} {params.filters} -uS 2>> {log} |\
            samtools sort -o {output.bam} >> {log} 2>&1
            samtools index {output.bam} >> {log} 2>&1
            """
    
    rule RemoveAdapters_p2:
        input: rules.RemoveAdapters_p1.output.bam
        output: f"{datadir + cln + noad}" + "{sample}.fastq"
        conda:
            f"{conda_envs}Clean.yaml"
        threads: config['threads']['AdapterRemoval']
        resources:
            mem_mb = low_memory_job
        params:
            script = srcdir('scripts/clipper.py')
        shell: 
            """
            python {params.script} --input {input} --output {output} --threads {threads}
            """

    rule QC_filter:
        input: rules.RemoveAdapters_p2.output
        output: 
            fq = f"{datadir + cln + qcfilt}" + "{sample}.fastq",
            html = f"{datadir + cln + qcfilt + html}" + "{sample}.fastp.html",
            json = f"{datadir + cln + qcfilt + json}" + "{sample}.fastp.json"
        conda: 
            f"{conda_envs}Clean.yaml"
        log:
            f"{logdir}" + "Cleanup_{sample}.log"
        benchmark:
            f"{logdir + bench}" + "Cleanup_{sample}.txt"
        threads: config['threads']['QC']
        resources:
            mem_mb = low_memory_job
        params:
            score = config['runparams']['qc_filter_iontorrent'],
            size = config['runparams']['qc_window_iontorrent'],
            length = config['runparams']['qc_min_readlength']
        shell:
            """
            fastp --thread {threads} -i {input} \
            --cut_right -M {params.score} -W {params.size} -l {params.length} \
            -o {output.fq} -h {output.html} -j {output.json} > {log} 2>&1
            """

if config["primer_file"] != "NONE":
    rule RemovePrimers:
        input: 
            fq = rules.QC_filter.output.fq,
            pr = rules.Prepare_ref_and_primers.output.prm,
            ref = rules.Prepare_ref_and_primers.output.ref
        output: 
            fq = f"{datadir + cln + prdir}" + "{sample}.fastq",
            ep = f"{datadir + prim}" + "{sample}_removedprimers.csv"
        conda: 
            f"{conda_envs}Clean.yaml"
        log:
            f"{logdir}" + "RemovePrimers_{sample}.log"
        benchmark:
            f"{logdir + bench}" + "RemovePrimers_{sample}.txt"
        threads: config['threads']['PrimerRemoval']
        resources:
            mem_mb = high_memory_job
        params:
            amplicontype = config["amplicon_type"]
        shell:
            """
            AmpliGone -i {input.fq} \
            -ref {input.ref} -pr {input.pr} \
            -o {output.fq} \
            -at {params.amplicontype} \
            --export-primers {output.ep} \
            -to \
            -t {threads}
            """
if config["primer_file"] == "NONE":
    rule RemovePrimers:
        input: rules.QC_filter.output.fq
        resources:
            mem_mb = 1
        output: 
            fq = f"{datadir + cln + prdir}" + "{sample}.fastq"
        shell:
            """
            cp {input} {output}
            """

rule QC_clean:
    input: rules.RemovePrimers.output.fq
    output: 
        html    =   f"{datadir + qc_post}" + "{sample}_fastqc.html",
        zip     =   f"{datadir + qc_post}" + "{sample}_fastqc.zip"
    conda: 
        f"{conda_envs}Clean.yaml"
    log:
        f"{logdir}" + "QC_clean_data_{sample}.log"
    benchmark:
        f"{logdir + bench}" + "QC_clean_data_{sample}.txt"
    threads: config['threads']['QC']
    resources:
        mem_mb = low_memory_job
    params:
        outdir = f"{datadir + qc_post}"
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
if config["platform"] == "illumina":
    rule Alignment:
        input:
            fq = rules.RemovePrimers.output.fq,
            ref = rules.Prepare_ref_and_primers.output.ref
        output:
            bam = f"{datadir + aln + bf}" + "{sample}.bam",
            index = f"{datadir + aln + bf}" + "{sample}.bam.bai"
        conda: 
            f"{conda_envs}Alignment.yaml"
        log:
            f"{logdir}" + "Alignment_{sample}.log"
        benchmark:
            f"{logdir + bench}" + "Alignment_{sample}.txt"
        threads: config['threads']['Alignments']
        resources:
            mem_mb = medium_memory_job
        params:
            mapthreads = config['threads']['Alignments'] - 1,
            filters = config["runparams"]["alignmentfilters"]
        shell: 
            """
            minimap2 -ax sr -t {params.mapthreads} {input.ref} {input.fq} 2>> {log} |\
            samtools view -@ {threads} {params.filters} -uS 2>> {log} |\
            samtools sort -o {output.bam} >> {log} 2>&1
            samtools index {output.bam} >> {log} 2>&1
            """

if config["platform"] == "nanopore":
    rule Alignment:
        input:
            fq = rules.RemovePrimers.output.fq,
            ref = rules.Prepare_ref_and_primers.output.ref
        output:
            bam = f"{datadir + aln + bf}" + "{sample}.bam",
            index = f"{datadir + aln + bf}" + "{sample}.bam.bai"
        conda: 
            f"{conda_envs}Alignment.yaml"
        log:
            f"{logdir}" + "Alignment_{sample}.log"
        benchmark:
            f"{logdir + bench}" + "Alignment_{sample}.txt"
        threads: config['threads']['Alignments']
        resources:
            mem_mb = medium_memory_job
        params:
            mapthreads = config['threads']['Alignments'] - 1,
            filters = config["runparams"]["alignmentfilters"]
        shell: 
            """
            minimap2 -ax map-ont -t {params.mapthreads} {input.ref} {input.fq} 2>> {log} |\
            samtools view -@ {threads} {params.filters} -uS 2>> {log} |\
            samtools sort -o {output.bam} >> {log} 2>&1
            samtools index {output.bam} >> {log} 2>&1
            """

if config["platform"] == "iontorrent":
    rule Alignment:
        input:
            fq = rules.RemovePrimers.output.fq,
            ref = rules.Prepare_ref_and_primers.output.ref
        output:
            bam = f"{datadir + aln + bf}" + "{sample}.bam",
            index = f"{datadir + aln + bf}" + "{sample}.bam.bai"
        conda: 
            f"{conda_envs}Alignment.yaml"
        log:
            f"{logdir}" + "Alignment_{sample}.log"
        benchmark:
            f"{logdir + bench}" + "Alignment_{sample}.txt"
        threads: config['threads']['Alignments']
        resources:
            mem_mb = medium_memory_job
        params:
            mapthreads = config['threads']['Alignments'] - 1,
            filters = config["runparams"]["alignmentfilters"]
        shell: 
            """
            minimap2 -ax sr -t {params.mapthreads} {input.ref} {input.fq} 2>> {log} |\
            samtools view -@ {threads} {params.filters} -uS 2>> {log} |\
            samtools sort -o {output.bam} >> {log} 2>&1
            samtools index {output.bam} >> {log} 2>&1
            """

rule Consensus:
    input: 
        bam = rules.Alignment.output.bam,
        gff = rules.Prepare_features_file.output.gff,
        ref = rules.Prepare_ref_and_primers.output.ref
    output: 
        cons_1 = f"{datadir + cons + seqs}" + "{sample}_cov_ge_1.fa",
        cons_5 = f"{datadir + cons + seqs}" + "{sample}_cov_ge_5.fa",
        cons_10 = f"{datadir + cons + seqs}" + "{sample}_cov_ge_10.fa",
        cons_50 = f"{datadir + cons + seqs}" + "{sample}_cov_ge_50.fa",
        cons_100 = f"{datadir + cons + seqs}" + "{sample}_cov_ge_100.fa",
        cov = f"{datadir + cons + covs}" + "{sample}_coverage.tsv",
        vcf_1 = f"{datadir + aln + vf}" + "{sample}_cov_ge_1.vcf",
        vcf_5 = f"{datadir + aln + vf}" + "{sample}_cov_ge_5.vcf",
        vcf_10 = f"{datadir + aln + vf}" + "{sample}_cov_ge_10.vcf",
        vcf_50 = f"{datadir + aln + vf}" + "{sample}_cov_ge_50.vcf",
        vcf_100 = f"{datadir + aln + vf}" + "{sample}_cov_ge_100.vcf",
        gff_1 = f"{datadir + cons + features}" + "{sample}_cov_ge_1.gff",
        gff_5 = f"{datadir + cons + features}" + "{sample}_cov_ge_5.gff",
        gff_10 = f"{datadir + cons + features}" + "{sample}_cov_ge_10.gff",
        gff_50 = f"{datadir + cons + features}" + "{sample}_cov_ge_50.gff",
        gff_100 = f"{datadir + cons + features}" + "{sample}_cov_ge_100.gff"
    params:
        mincov = "1 5 10 50 100",
        outdir = f"{datadir + cons + seqs}",
        vcfdir = f"{datadir + aln + vf}",
        gffdir = f"{datadir + cons + features}"
    conda: 
        f"{conda_envs}Consensus.yaml"
    log:
        f"{logdir}" + "Consensus_{sample}.log"
    benchmark:
        f"{logdir + bench}" + "Consensus_{sample}.txt"
    threads: config['threads']['Consensus']
    resources:
        mem_mb = medium_memory_job
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

rule Concat_Seqs:
    input: 
        cov1 =  expand( "{p}{sample}_cov_ge_1.fa",
                    p = f"{datadir + cons + seqs}",
                    sample = SAMPLES),
        cov5 =  expand( "{p}{sample}_cov_ge_5.fa",
                    p = f"{datadir + cons + seqs}",
                    sample = SAMPLES),
        cov10 = expand( "{p}{sample}_cov_ge_10.fa",
                    p = f"{datadir + cons + seqs}",
                    sample = SAMPLES),
        cov50 = expand( "{p}{sample}_cov_ge_50.fa",
                    p = f"{datadir + cons + seqs}",
                    sample = SAMPLES),
        cov100 =expand( "{p}{sample}_cov_ge_100.fa",
                    p = f"{datadir + cons + seqs}",
                    sample = SAMPLES),
    output: 
        cov1 = f"{res + seqs}" + "concat_cov_ge_1.fasta",
        cov5 = f"{res + seqs}" + "concat_cov_ge_5.fasta",
        cov10 = f"{res + seqs}" + "concat_cov_ge_10.fasta",
        cov50 = f"{res + seqs}" + "concat_cov_ge_50.fasta",
        cov100 = f"{res + seqs}" + "concat_cov_ge_100.fasta",
    threads: 1
    resources:
        mem_mb = low_memory_job
    shell:
        """
        cat {input.cov1} >> {output.cov1}
        cat {input.cov5} >> {output.cov5}
        cat {input.cov10} >> {output.cov10}
        cat {input.cov50} >> {output.cov50}
        cat {input.cov100} >> {output.cov100}
        """


rule VCF_to_TSV:
    input: 
        cov1 = f"{datadir + aln + vf}" + "{sample}_cov_ge_1.vcf",
        cov5 = f"{datadir + aln + vf}" + "{sample}_cov_ge_5.vcf",
        cov10 = f"{datadir + aln + vf}" + "{sample}_cov_ge_10.vcf",
        cov50 = f"{datadir + aln + vf}" + "{sample}_cov_ge_50.vcf",
        cov100 = f"{datadir + aln + vf}" + "{sample}_cov_ge_100.vcf",
    output: 
        cov1 = temp(f"{datadir + aln + vf}" + "{sample}_cov_ge_1.tsv"),
        cov5 = temp(f"{datadir + aln + vf}" + "{sample}_cov_ge_5.tsv"),
        cov10 = temp(f"{datadir + aln + vf}" + "{sample}_cov_ge_10.tsv"),
        cov50 = temp(f"{datadir + aln + vf}" + "{sample}_cov_ge_50.tsv"),
        cov100 = temp(f"{datadir + aln + vf}" + "{sample}_cov_ge_100.tsv"),
    conda: f"{conda_envs}Mutations.yaml"
    threads: config['threads']['Index']
    resources:
        mem_mb = low_memory_job
    log:
        f"{logdir}" + "vcf_to_tsv_{sample}.log"
    shell:
        """
        bcftools query {input.cov1} -f '{wildcards.sample}\t%CHROM\t%POS\t%REF\t%ALT\t%DP\n' -e 'ALT="N"' > {output.cov1} 2>> {log}
        bcftools query {input.cov5} -f '{wildcards.sample}\t%CHROM\t%POS\t%REF\t%ALT\t%DP\n' -e 'ALT="N"' > {output.cov5} 2>> {log}
        bcftools query {input.cov10} -f '{wildcards.sample}\t%CHROM\t%POS\t%REF\t%ALT\t%DP\n' -e 'ALT="N"' > {output.cov10} 2>> {log}
        bcftools query {input.cov50} -f '{wildcards.sample}\t%CHROM\t%POS\t%REF\t%ALT\t%DP\n' -e 'ALT="N"' > {output.cov50} 2>> {log}
        bcftools query {input.cov100} -f '{wildcards.sample}\t%CHROM\t%POS\t%REF\t%ALT\t%DP\n' -e 'ALT="N"' > {output.cov100} 2>> {log}
        """

rule Concat_TSV_coverages:
    input:
        cov1 = expand(
            "{p}{sample}_cov_ge_1.tsv",
            p = datadir + aln + vf,
            sample = SAMPLES
            ),
        cov5 = expand(
            "{p}{sample}_cov_ge_5.tsv",
            p = datadir + aln + vf,
            sample = SAMPLES
            ),
        cov10 = expand(
            "{p}{sample}_cov_ge_10.tsv",
            p = datadir + aln + vf,
            sample = SAMPLES
            ),
        cov50 = expand(
            "{p}{sample}_cov_ge_50.tsv",
            p = datadir + aln + vf,
            sample = SAMPLES
            ),
        cov100 = expand(
            "{p}{sample}_cov_ge_100.tsv",
            p = datadir + aln + vf,
            sample = SAMPLES
            )
    output:
        cov1 = f"{res + muts}concat_mutations_cov_ge_1.tsv",
        cov5 = f"{res + muts}concat_mutations_cov_ge_5.tsv",
        cov10 = f"{res + muts}concat_mutations_cov_ge_10.tsv",
        cov50 = f"{res + muts}concat_mutations_cov_ge_50.tsv",
        cov100 = f"{res + muts}concat_mutations_cov_ge_100.tsv",
    log:
        f"{logdir}concat_tsv.log"
    threads: 1
    resources:
        mem_mb = low_memory_job
    run:
        shell("echo -e 'Sample\tReference_Chromosome\tPosition\tReference\tAlternative\tDepth' > {output.cov1} 2> {log}")
        shell("echo -e 'Sample\tReference_Chromosome\tPosition\tReference\tAlternative\tDepth' > {output.cov5} 2> {log}")
        shell("echo -e 'Sample\tReference_Chromosome\tPosition\tReference\tAlternative\tDepth' > {output.cov10} 2> {log}")
        shell("echo -e 'Sample\tReference_Chromosome\tPosition\tReference\tAlternative\tDepth' > {output.cov50} 2> {log}")
        shell("echo -e 'Sample\tReference_Chromosome\tPosition\tReference\tAlternative\tDepth' > {output.cov100} 2> {log}")
        shell("cat {input.cov1} >> {output.cov1} 2>> {log}")
        shell("cat {input.cov5} >> {output.cov5} 2>> {log}")
        shell("cat {input.cov10} >> {output.cov10} 2>> {log}")
        shell("cat {input.cov50} >> {output.cov50} 2>> {log}")
        shell("cat {input.cov100} >> {output.cov100} 2>> {log}")


rule Get_Breadth_of_coverage:
    input:
        reference = rules.Prepare_ref_and_primers.output.ref,
        coverage = rules.Consensus.output.cov,
    output:
        temp(f"{datadir + boc}" + "{sample}.tsv")
    conda:
        f"{conda_envs}Consensus.yaml"
    threads: 1
    resources:
        mem_mb = low_memory_job
    params:
        script = srcdir("scripts/boc.py")
    shell:
        """
        python {params.script} {input.reference} {wildcards.sample} {input.coverage} {output}
        """

rule concat_boc:
    input:
        pct = expand("{p}{sample}.tsv",
                    p = f"{datadir + boc}",
                    sample = SAMPLES)
    output: f"{res}Width_of_coverage.tsv"
    threads: 1
    resources:
        mem_mb = low_memory_job
    shell:
        """
        echo -e "Sample_name\tWidth_at_mincov_1\tWidth_at_mincov_5\tWidth_at_mincov_10\tWidth_at_mincov_50\tWidth_at_mincov_100" > {output}
        cat {input} >> {output}
        """
if config["primer_file"] != "NONE":
    rule calculate_amplicon_cov:
        input: 
            pr = rules.RemovePrimers.output.ep,
            cov = rules.Consensus.output.cov
        output: 
            ampcov = f"{datadir + prim}" + "{sample}_ampliconcoverage.csv"
        threads: 1
        resources:
            mem_mb = low_memory_job
        conda:
            f"{conda_envs}Clean.yaml"
        params:
            script = srcdir("scripts/amplicon_covs.py")
        shell:
            """
            python {params.script} \
            --primers {input.pr} \
            --coverages {input.cov} \
            --key {wildcards.sample} \
            --output {output.ampcov}
            """
    
    rule concat_amplicon_cov:
        input: expand(f"{datadir + prim}" + "{sample}_ampliconcoverage.csv", sample = SAMPLES)
        output: f"{res}Amplicon_coverage.csv"
        threads: 1
        resources:
            mem_mb = low_memory_job
        conda:
            f"{conda_envs}Clean.yaml"
        params:
            script = srcdir("scripts/concat_amplicon_covs.py")
        shell:
            """
            python {params.script} --output {output} --input {input}
            """
    
if config['platform'] == "illumina":
    rule MultiQC_report:
        input: 
            expand("{p}{sample}_{read}_fastqc.zip",
                p = f"{datadir + qc_pre}",
                sample = SAMPLES,
                read = "R1 R2".split()
                ),
            expand("{p}{sample}_fastqc.zip",
                p = f"{datadir + qc_post}",
                sample = SAMPLES
                )
        output: 
            f"{res}multiqc.html",
            expand("{p}multiqc_{program}.txt",
                p = f"{res + mqc_data}",
                program = ['fastqc'])
        conda:
            f"{conda_envs}Clean.yaml"
        log:
            f"{logdir}" + "MultiQC_report.log"
        benchmark:
            f"{logdir + bench}" + "MultiQC_report.txt"
        threads: 1
        resources:
            mem_mb = medium_memory_job
        params:
            conffile = srcdir('files/multiqc_config.yaml'),
            outdir = f"{res}"
        shell:
            """
            multiqc -d --force --config {params.conffile} -o {params.outdir} -n multiqc.html {input} > {log} 2>&1
            """

if config['platform'] == "nanopore" or config['platform'] == "iontorrent":
    rule MultiQC_report:
        input: 
            expand("{p}{sample}_fastqc.zip",
                p = f"{datadir + qc_pre}",
                sample = SAMPLES
                ),
            expand("{p}{sample}_fastqc.zip",
                p = f"{datadir + qc_post}",
                sample = SAMPLES
                )
        output: 
            f"{res}multiqc.html",
            expand("{p}multiqc_{program}.txt",
                p = f"{res + mqc_data}",
                program = ['fastqc'])
        conda:
            f"{conda_envs}Clean.yaml"
        log:
            f"{logdir}" + "MultiQC_report.log"
        benchmark:
            f"{logdir + bench}" + "MultiQC_report.txt"
        threads: 1
        resources:
            mem_mb = medium_memory_job
        params:
            conffile = srcdir('files/multiqc_config.yaml'),
            outdir = f"{res}"
        shell:
            """
            multiqc -d --force --config {params.conffile} -o {params.outdir} -n multiqc.html {input} > {log} 2>&1
            """