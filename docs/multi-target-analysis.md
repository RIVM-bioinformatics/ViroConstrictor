# Multi-Target Analysis

Multi-target analysis allows you to apply different settings to each sample within a single analysis run. This is useful when analyzing samples with different viral targets, experimental conditions, or analysis parameters.

## When to Use Multi-Target Analysis

Use multi-target analysis when you have:

* Multiple viral targets in the same batch
* Different primer sets for different samples
* Varying analysis parameters (coverage thresholds, mismatch rates)
* Mixed experimental protocols

## Samplesheet Format

Multi-target analysis requires a samplesheet in Excel (`.xlsx`), CSV (`.csv`), or TSV (`.tsv`) format. You can download an example Excel spreadsheet [here](samples.xlsx).

### Required Columns

!!! attention "Mandatory Columns"
    The following columns are required and analysis will fail without them:

    * **Sample**: Must match FASTQ filenames exactly
    * **Virus**: Viral target name (triggers presets if available)
    * **Reference**: Path to reference FASTA file

### Optional Columns

All other columns are optional and will use default values if not specified:

| Column | Purpose | Example Values |
|--------|---------|----------------|
| `Primers` | Primer file path | `/path/to/primers.fasta` or `NONE` |
| `Features` | GFF file path | `/path/to/features.gff` or `NONE` |
| `Match-ref` | Enable reference matching | `TRUE` or `FALSE` |
| `Segmented` | Segmented virus analysis | `TRUE` or `FALSE` |
| `min-coverage` | Minimum coverage threshold | `30`, `50`, etc. |
| `primer-mismatch-rate` | Primer mismatch tolerance | `0.1`, `0.15`, etc. |

<!-- TODO: Add guidance on column naming conventions and case sensitivity -->

### Example Samplesheet

| Sample    | Virus       | Match-ref | Segmented | Primers                              | Reference                                 | Features           | min-coverage | primer-mismatch-rate |
|-----------|-------------|-----------|-----------|--------------------------------------|-------------------------------------------|--------------------|--------------|----------------------|
| Sample_1  | SARS-CoV-2  | FALSE     | FALSE     | /path/to/sars_cov_2_primers.fasta    | /path/to/sars_cov_2_reference.fasta       | /path/to/sars_cov_2.gff | 50           | 0.1                  |
| Sample_2  | Measles     | FALSE     | FALSE     | /path/to/Measles_primers.fasta        | /path/to/measles_reference.fasta          | /path/to/measles.gff | 30           | 0.15                 |
| Sample_3  | Influenza_A | TRUE      | TRUE      | /path/to/Influenza-A_primers.fasta    | /path/to/Influenza-A_reference.fasta       | NONE               | 30           | 0.1                  |
| Sample_4  | HPV         | TRUE      | FALSE     | NONE                                 | /path/to/HPV_reference.fasta              | NONE               | 10           | 0.1                  |

## Important Notes

### Sample Naming
!!! warning "Exact Filename Matching Required"
    The `Sample` column must correspond exactly to your FASTQ filenames in the input directory. ViroConstrictor will alert you if names don't match to prevent unexpected results.

### Viral Target Names
- The `Virus` column serves as the viral target name for analysis
- Names trigger analysis presets if ViroConstrictor recognizes them
- If no preset exists, default settings are applied
- See [presets documentation](presets.md) for available presets
- Use `--disable-presets` flag to force default settings for all samples

<!-- TODO: Add examples of common sample naming issues and solutions -->

## Running Multi-Target Analysis

### Basic Command

!!! example "Multi-Target Analysis Command"
    ```bash
    viroconstrictor \
        --input /path/to/FASTQ-files \
        --output /path/to/output-folder \
        --samplesheet /path/to/samplesheet.xlsx \
        --platform nanopore \
        --amplicon-type end-to-end
    ```

### Combining Samplesheet with Command-Line Parameters

You can combine samplesheet settings with run-wide command-line parameters. This is useful when most parameters are the same across samples, but a few need to be customized per sample.

!!! example "Mixed Parameter Approach"
    ```bash
    # Common platform and amplicon type for all samples
    # Individual coverage thresholds specified in samplesheet
    viroconstrictor \
        --input /path/to/FASTQ-files \
        --output /path/to/output-folder \
        --samplesheet /path/to/samplesheet.xlsx \
        --platform illumina \
        --amplicon-type end-to-end \
        --primer-mismatch-rate 0.15
    ```

**Parameter priority**: Samplesheet values override command-line parameters for individual samples.

<!-- TODO: Clarify exactly which parameters can be mixed and which cannot -->
