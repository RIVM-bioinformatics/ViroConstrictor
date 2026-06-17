rule align_before_trueconsense:
    input:
        fq=f"{datadir}{wc_folder}{cln}{prdir}" "{sample}.fastq",  #rules.ampligone.output.fq,
        ref=rules.prepare_refs.output,
    output:
        bam=f"{datadir}{wc_folder}{aln}{bf}" "{sample}.bam",
        index=f"{datadir}{wc_folder}{aln}{bf}" "{sample}.bam.bai",
    conda:
        workflow_environment_path("Alignment.yaml")
    container:
        f"{container_base_path}/viroconstrictor_alignment_{get_hash('Alignment')}.sif"
    log:
        f"{logdir}{get_rule_name()}/Alignment_" "{Virus}.{RefID}.{sample}.log",
    threads: config["threads"]["Alignments"]
    resources:
        mem_mb=medium_memory_job,
        runtime=medium_runtime_job,
    params:
        mapthreads=config["threads"]["Alignments"] - 1,
        align_bin = lambda wc: get_binary(
            stage_identifier=VC_STAGE,
            preset_name=SAMPLES[wc.sample]["PRESET"],
            task="alignment",
            platform=config["platform"],
        ),
        align_flags = lambda wc: get_flags(
            stage_identifier=VC_STAGE,
            preset_name=SAMPLES[wc.sample]["PRESET"],
            task="alignment",
            platform=config["platform"],
        ),
        filter_bin = lambda wc: get_binary(
            stage_identifier=VC_STAGE,
            preset_name=SAMPLES[wc.sample]["PRESET"],
            task="alignment_filtering",
            platform=config["platform"],
        ),
        filter_flags = lambda wc: get_flags(
            stage_identifier=VC_STAGE,
            preset_name=SAMPLES[wc.sample]["PRESET"],
            task="alignment_filtering",
            platform=config["platform"],
        ),
    shell:
        """
        {params.align_bin} {params.align_flags} -t {params.mapthreads} {input.ref} {input.fq} 2>> {log} |\
        {params.filter_bin} view -@ {threads} {params.filter_flags} -uS 2>> {log} |\
        {params.filter_bin} sort -o {output.bam} >> {log} 2>&1
        {params.filter_bin} index {output.bam} >> {log} 2>&1
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
        workflow_environment_path("Consensus.yaml")
    container:
        f"{container_base_path}/viroconstrictor_consensus_{get_hash('Consensus')}.sif"
    log:
        f"{logdir}{get_rule_name()}/Consensus_" "{Virus}.{RefID}.{sample}.log",
    threads: config["threads"]["Consensus"]
    resources:
        mem_mb=medium_memory_job,
        runtime=high_runtime_job,
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


rule Translate_AminoAcids:
    input:
        seq=rules.trueconsense.output.cons,
        gff=rules.trueconsense.output.gff,
    output:
        f"{datadir}{wc_folder}{amino}" "{sample}/aa.faa",
    conda:
        workflow_environment_path("ORF_analysis.yaml")
    container:
        f"{container_base_path}/viroconstrictor_orf_analysis_{get_hash('ORF_analysis')}.sif"
    resources:
        mem_mb=low_memory_job,
        runtime=low_runtime_job,
    log:
        f"{logdir}{get_rule_name()}/Translate_AA_" "{Virus}.{RefID}.{sample}.log",
    params:
        binary = lambda wc: get_binary(
            stage_identifier=VC_STAGE,
            preset_name=SAMPLES[wc.sample]["PRESET"],
            task="amino_acid_extraction",
            platform=config["platform"],
        ),
        flags = lambda wc: get_flags(
            stage_identifier=VC_STAGE,
            preset_name=SAMPLES[wc.sample]["PRESET"],
            task="amino_acid_extraction",
            platform=config["platform"],
        ),
    shell:
        """
        {params.binary} -i {input.seq} -gff {input.gff} -o {output} -n {wildcards.sample} {params.flags} >> {log} 2>&1
        """
