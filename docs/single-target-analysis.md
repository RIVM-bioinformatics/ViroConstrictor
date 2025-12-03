# Single-Target Analysis

Single-target analysis applies the same settings to all samples in your dataset. This is the simplest way to run ViroConstrictor when all your samples share the same viral target and experimental conditions.

## Required Parameters

To run a single-target analysis, you need to provide the following inputs:

### Input/Output Directories
- **Input directory**: `--input` / `-i` - Folder containing your FASTQ files
- **Output directory**: `--output` / `-o` - Desired location for results and data

### Reference and Annotations
- **Reference sequence**: `--reference` / `-ref` - Path to reference FASTA file
- **Primers**: `--primers` / `-pr` - Primer sequences in FASTA or BED format
- **Genomic features**: `--features` / `-gff` - GFF3 file with genomic annotations

<!-- TODO: Add guidance on what to do when reference files are missing or unavailable -->

!!! info "No Primers Used?"
    If your sequencing protocol didn't use primers, specify `--primers NONE`

!!! info "No GFF File Available?"
    If you don't have a matching GFF file, use `--features NONE`. ViroConstrictor will attempt to infer genomic features, but amino acid translations won't be included in results. For amino acid translations, provide a matching GFF file.

### Sequencing Configuration
- **Platform**: `--platform` - Sequencing technology used (`nanopore`, `illumina`, or `iontorrent`)
- **Amplicon type**: `--amplicon-type` / `-at` - Type of amplicons in your data (see [amplicon documentation](amplicons.md))

!!! info "Primer-Free Protocols"
    If you specified `--primers NONE`, the amplicon type setting will be ignored but still needs to be provided.

### Analysis Settings  
- **Viral target**: `--target` / `--preset` - Name of the viral target being analyzed

!!! attention "Analysis Presets"
    The viral target name will trigger analysis presets if ViroConstrictor recognizes it. If no preset exists, default settings are used. See [presets documentation](presets.md) for details. Use `--disable-presets` to force default settings.

## Optional Parameters

You can fine-tune the analysis with these optional settings:

- **Minimum coverage**: `--min-coverage` / `-mc` - Minimum coverage threshold (default: 30)
- **Primer mismatch rate**: `--primer-mismatch-rate` / `-pmr` - Maximum allowed primer mismatches (default: 0.1 = 10%)

<!-- TODO: Explain the impact of changing these parameters on analysis results -->

## Example Command

!!! example "Basic Single-Target Analysis"
    ```bash
    viroconstrictor \
        --input /path/to/FASTQ-files \
        --output /path/to/output-folder \
        --reference /path/to/reference.fasta \
        --primers /path/to/primers.fasta \
        --features /path/to/features.gff \
        --platform nanopore \
        --amplicon-type end-to-end \
        --target SARS-CoV-2
    ```

!!! example "Analysis with Custom Parameters"
    ```bash
    viroconstrictor \
        --input /path/to/FASTQ-files \
        --output /path/to/output-folder \
        --reference /path/to/reference.fasta \
        --primers /path/to/primers.fasta \
        --features /path/to/features.gff \
        --platform illumina \
        --amplicon-type end-to-end \
        --target HPV16 \
        --min-coverage 50 \
        --primer-mismatch-rate 0.15
    ```

<!-- TODO: Add example for primer-free analysis -->

