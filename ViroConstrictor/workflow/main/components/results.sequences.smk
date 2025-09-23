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
        workflow_environment_path("Alignment.yaml")
    container:
        f"{container_base_path}/viroconstrictor_alignment_{get_hash('Alignment')}.sif"
    log:
        f"{logdir}Alignment_" "{Virus}.{RefID}.{sample}.log",
    benchmark:
        f"{logdir}{bench}Alignment_" "{Virus}.{RefID}.{sample}.txt"
    threads: config["threads"]["Alignments"]
    resources:
        mem_mb=medium_memory_job,
        runtime=55,
    params:
        mapthreads=config["threads"]["Alignments"] - 1,
        mm2_alignment_preset=get_alignment_flags,
        minimap2_base_setting=lambda wc: get_preset_parameter(
            preset_name=SAMPLES[wc.sample]["PRESET"],
            parameter_name="Minimap2_Settings_Base",
        ),
        minimap2_extra_setting=lambda wc: get_preset_parameter(
            preset_name=SAMPLES[wc.sample]["PRESET"],
            parameter_name="Minimap2_Settings",
        ),
        minimap2_alignmentparams=lambda wc: get_preset_parameter(
            preset_name=SAMPLES[wc.sample]["PRESET"],
            parameter_name=f"Minimap2_AlignmentParams_{config['platform']}",
        ),
        samtools_standard_filters=lambda wc: get_preset_parameter(
            preset_name=SAMPLES[wc.sample]["PRESET"],
            parameter_name="Samtools_Filters_Base",
        ),
        samtools_extra_filters=lambda wc: get_preset_parameter(
            preset_name=SAMPLES[wc.sample]["PRESET"],
            parameter_name=f"Samtools_Filters_{config['platform']}",
        ),
    shell:
        """
        minimap2 {params.mm2_alignment_preset} {params.minimap2_base_setting} {params.minimap2_extra_setting} {params.minimap2_alignmentparams} -t {params.mapthreads} {input.ref} {input.fq} 2>> {log} |\
        samtools view -@ {threads} {params.samtools_standard_filters} {params.samtools_extra_filters} -uS 2>> {log} |\
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
        workflow_environment_path("Consensus.yaml")
    container:
        f"{container_base_path}/viroconstrictor_consensus_{get_hash('Consensus')}.sif"
    log:
        f"{logdir}Consensus_" "{Virus}.{RefID}.{sample}.log",
    benchmark:
        f"{logdir}{bench}Trueconsense_" "{Virus}.{RefID}.{sample}.txt"
    threads: config["threads"]["Consensus"]
    resources:
        mem_mb=medium_memory_job,
        runtime=55,
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
        runtime=55,
    log:
        f"{logdir}Translate_AA_" "{Virus}.{RefID}.{sample}.log",
    benchmark:
        f"{logdir}{bench}Translate_AA_" "{Virus}.{RefID}.{sample}.txt"
    params:
        feature_type=lambda wc: get_preset_parameter(
            preset_name=SAMPLES[wc.sample]["PRESET"],
            parameter_name="AminoExtract_FeatureType",
        ),
        aminoextract_settings=lambda wc: get_preset_parameter(
            preset_name=SAMPLES[wc.sample]["PRESET"],
            parameter_name=f"AminoExtract_Settings",
        ),
    shell:
        """
        AminoExtract -i {input.seq} -gff {input.gff} -o {output} -ft {params.feature_type} -n {wildcards.sample} {params.aminoextract_settings} >> {log} 2>&1
        """
