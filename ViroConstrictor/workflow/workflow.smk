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
    SAMPLES = yaml.safe_load(sample_sheet_file)

primerfile = config["primer_file"]

if primerfile == "NONE":
    primerfile = srcdir("files/empty.primers")
primers_extension = primerfile.split('.')[-1]

reffile = config["reference_file"]
ref_basename = os.path.splitext(os.path.basename(reffile))[0]

features_file = config["features_file"]

mincov = 30

def construct_all_rule(_wildcards):
    files = []
    files.append(f"{res}multiqc.html")
    files.append(f"{res}consensus.fasta")
    files.append(f"{res}mutations.tsv")
    files.append(f"{res}Width_of_coverage.tsv")

    if config["primer_file"] == "NONE":
        files.append(f"{res}Amplicon_coverage.csv")

    return files

rule all:
    input: construct_all_rule


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

rule Prepare_ref:
    input: ref = reffile,
    output:
        ref = f"{datadir}{refdir}{ref_basename}.fasta",
        refindex = f"{datadir}{refdir}{ref_basename}.fasta.fai",
    conda: f"{conda_envs}Alignment.yaml"
    threads: config['threads']['Index']
    resources: mem_mb = low_memory_job
    shell:
        """
        cat {input.ref} | seqkit replace -p "\-" -s -r "N" > {output.ref}
        samtools faidx {output.ref} -o {output.refindex}
            """

rule CopyPrimersToData:
    input: primerfile
    output: f"{datadir}{prim}primers.{primers_extension}"
    threads: 1
    resources: mem_mb = low_memory_job
    shell: "cp {input} {output}"

rule PrimersToBED:
    input:
        prm = rules.CopyPrimersToData.output,
        ref = rules.Prepare_ref.output.ref
    output: primer_bed = f"{datadir}{prim}""primers.bed"
    conda: f"{conda_envs}Clean.yaml"
    params: pr_mm_rate = config["primer_mismatch_rate"]
    log: f"{logdir}PrimersToBED.log"
    threads: 1
    shell:
        """
        python -m AmpliGone.fasta2bed \
            --primers {input.prm} \
            --reference {input.ref} \
            --output {output.primer_bed} \
            --primer-mismatch-rate {params.pr_mm_rate} > {log}
        """

if config["features_file"] != "NONE":
    rule Prepare_features_file:
        input: features_file
        output: gff = temp(f"{datadir}{features}reference_features.gff")
        threads: config['threads']['Index']
        resources: mem_mb = low_memory_job
        shell:
            """
            cat {input} >> {output}
            """
else:
    rule Prepare_features_file:
        input: rules.Prepare_ref.output.ref
        output:
            AA  = temp(f"{datadir}{refdir}{ref_basename}_ORF_AA.fa"),
            NT  = temp(f"{datadir}{refdir}{ref_basename}_ORF_NT.fa"),
            gff = temp(f"{datadir}{refdir}{ref_basename}_annotation.gff"),
        conda: f"{conda_envs}ORF_analysis.yaml"
        threads: config['threads']['Index']
        resources: mem_mb = medium_memory_job
        log: f"{logdir}ORF_Analysis.log"
        params:
            procedure = "meta",
            output_format = "gff",
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

    rule RemoveAdapters_p1:
        input:
            ref = rules.Prepare_ref.output.ref,
            fq  = lambda wildcards: (SAMPLES[wildcards.sample][i]
                                for i in ("R1", "R2")
                                )
        output:
            bam     = f"{datadir}{cln}{raln}""{sample}.bam",
            index   = f"{datadir}{cln}{raln}""{sample}.bam.bai"
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

    rule RemoveAdapters_p2:
        input: rules.RemoveAdapters_p1.output.bam
        output: f"{datadir}{cln}{noad}""{sample}.fastq"
        conda:
            f"{conda_envs}Clean.yaml"
        threads: config['threads']['AdapterRemoval']
        resources: mem_mb = low_memory_job
        params: script = srcdir('scripts/clipper.py')
        shell:
            """
            python {params.script} --input {input} --output {output} --threads {threads}
            """

    rule QC_filter:
        input: rules.RemoveAdapters_p2.output
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
            score = config['runparams']['qc_filter_illumina'],
            size = config['runparams']['qc_window_illumina'],
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

if config["platform"] == "nanopore":
    rule QC_raw:
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
            output_dir  =   f"{datadir}{qc_pre}",
            script = srcdir("scripts/fastqc_wrapper.sh")
        shell:
            """
            bash {params.script} {input} {params.output_dir} {output.html} {output.zip} {log}
            """

    rule RemoveAdapters_p1:
        input:
            ref = rules.Prepare_ref.output.ref,
            fq  = lambda wildcards: SAMPLES[wildcards.sample]
        output:
            bam     = f"{datadir}{cln}{raln}""{sample}.bam",
            index   = f"{datadir}{cln}{raln}""{sample}.bam.bai"
        conda: f"{conda_envs}Alignment.yaml"
        log: f"{logdir}""RemoveAdapters_p1_{sample}.log"
        benchmark: f"{logdir}{bench}"+ "RemoveAdapters_p1_{sample}.txt"
        threads: config['threads']['Alignments']
        resources: mem_mb = medium_memory_job
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
        output: f"{datadir}{cln}{noad}""{sample}.fastq"
        conda: f"{conda_envs}Clean.yaml"
        threads: config['threads']['AdapterRemoval']
        resources: mem_mb = low_memory_job
        params: script = srcdir('scripts/clipper.py')
        shell:
            """
            python {params.script} --input {input} --output {output} --threads {threads}
            """

    rule QC_filter:
        input: rules.RemoveAdapters_p2.output
        output:
            fq   = f"{datadir}{cln}{qcfilt}""{sample}.fastq",
            html = f"{datadir}{cln}{qcfilt}{html}""{sample}.fastp.html",
            json = f"{datadir}{cln}{qcfilt}{json}""{sample}.fastp.json"
        conda: f"{conda_envs}Clean.yaml"
        log: f"{logdir}""Cleanup_{sample}.log"
        benchmark: f"{logdir}{bench}""Cleanup_{sample}.txt"
        threads: config['threads']['QC']
        resources: mem_mb = low_memory_job
        params:
            score = config['runparams']['qc_filter_nanopore'],
            size = config['runparams']['qc_window_nanopore'],
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

if config["platform"] == "iontorrent":
    rule QC_raw:
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

    rule RemoveAdapters_p1:
        input:
            ref = rules.Prepare_ref.output.ref,
            fq = lambda wildcards: SAMPLES[wildcards.sample]
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
        output: f"{datadir}{cln}{noad}""{sample}.fastq"
        conda: f"{conda_envs}Clean.yaml"
        threads: config['threads']['AdapterRemoval']
        resources: mem_mb = low_memory_job
        params: script = srcdir('scripts/clipper.py')
        shell:
            """
            python {params.script} --input {input} --output {output} --threads {threads}
            """

    rule QC_filter:
        input: rules.RemoveAdapters_p2.output
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
            score = config['runparams']['qc_filter_iontorrent'],
            size = config['runparams']['qc_window_iontorrent'],
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

ruleorder: AmpliGone > MoveFastq

rule AmpliGone:
    input:
        fq = rules.QC_filter.output.fq,
        pr = f"{datadir}{prim}primers.bed",
        ref = rules.Prepare_ref.output.ref
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

rule MoveFastq:
    input: rules.QC_filter.output.fq
    output: f"{datadir}{cln}{prdir}""{sample}.fastq"
    threads: 1
    resources: mem_mb = low_memory_job
    shell:
        """
        cp {input} {output}
        """

rule QC_clean:
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
    if config["platform"] in ["illumina", "iontorrent"]:
        return "-ax sr"
    elif config["platform"] == "nanopore":
        return "-ax map-ont -E2,0 -O8,24 -A4 -B4"

rule Alignment:
    input:
        fq = f"{datadir}{cln}{prdir}""{sample}.fastq",
        ref = rules.Prepare_ref.output.ref
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

rule Consensus:
    input:
        bam = rules.Alignment.output.bam,
        gff = rules.Prepare_features_file.output.gff,
        ref = rules.Prepare_ref.output.ref
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

rule Concat_Seqs:
    input:
        expand(
            f"{datadir}{cons}{seqs}""{sample}_cov_ge_"f"{mincov}.fa",
            sample = SAMPLES,
            )
    output: f"{res}consensus.fasta",
    threads: 1
    resources: mem_mb = low_memory_job
    shell: "cat {input} >> {output}"


rule VCF_to_TSV:
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

rule Concat_TSV_coverages:
    input: tsvs = expand(f"{datadir}{aln}{vf}""{sample}.tsv", sample = SAMPLES),
    output: tsv = f"{res}mutations.tsv",
    log: f"{logdir}concat_tsv.log"
    threads: 1
    resources: mem_mb = low_memory_job
    run:
        shell("echo -e 'Sample\tReference_Chromosome\tPosition\tReference\tAlternative\tDepth' > {output.tsv} 2> {log}")
        shell("cat {input.tsvs} >> {output.tsv} 2>> {log}")

rule Get_Breadth_of_coverage:
    input:
        reference = rules.Prepare_ref.output.ref,
        coverage = rules.Consensus.output.cov,
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
        pr = rules.AmpliGone.output.ep,
        cov = rules.Consensus.output.cov
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

def construct_MultiQC_input(_wildcards):
    if config['platform'] == "nanopore" or config['platform'] == "iontorrent":
        pre = expand(
            f"{datadir}{qc_pre}""{sample}_fastqc.zip",
            sample = SAMPLES,
            )
    elif config['platform'] == "illumina":
        pre = expand(
            f"{datadir}{qc_pre}""{sample}_{read}_fastqc.zip",
            sample = SAMPLES,
            read = "R1 R2".split()
            )
    else:
        raise ValueError(f"Platform {config['platform']} not recognised. Choose one of [illumina, nanopore, iontorrent].")
    post = expand(f"{datadir}{qc_post}""{sample}_fastqc.zip", sample = SAMPLES)
    return pre + post

rule MultiQC_report:
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
