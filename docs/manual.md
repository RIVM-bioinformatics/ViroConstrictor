# Starting an Analysis

Once ViroConstrictor is [installed](installation.md), you can use it to analyze your NGS FASTQ data.

You can run an analysis in two ways:

1. A [single-target](#run-a-single-target-analysis) analysis with run-wide parameters (the same settings are applied to every sample in the analysis).
2. A [multi-target](#run-a-multi-target-analysis) analysis where you provide a samplesheet (in Excel, CSV, or TSV format) that contains information for running the analysis with different settings for each sample.

Please see `viroconstrictor -h` for all command-line options.

---

## Overview of Command-Line Options

Below is a brief summary of the available command-line options and their meanings.

| Command                              | Argument                          | Explanation |
|--------------------------------------|-----------------------------------|-------------|
| `--input` /<br>`-i`                   | [input directory]                 | The folder containing your input FASTQ files. Both absolute and relative paths are supported. |
| `--output` /<br>`-o`                  | [output directory]                | The desired output folder for data and results. If the folder does not exist, it will be created. |
| `--reference` /<br>`-ref`             | [reference FASTA]                 | The input reference genome sequence in FASTA format. |
| `--primers` /<br>`-pr`                | [primer FASTA] / `NONE`           | The primer sequences in FASTA format.<br>Use `NONE` if your sequencing protocol does not use primers. |
| `--features` /<br>`-gff`              | [GFF file] / `NONE`               | The GFF3 file containing genomic features (e.g., open reading frames) that match the given reference FASTA.<br>Use `NONE` if you do not have a GFF3 file. |
| `--platform`                         | `nanopore` / `illumina` / `iontorrent` | The sequencing platform used to generate the dataset. The default is `nanopore`. |
| `--amplicon-type` /<br>`-at`          | `end-to-end` / `end-to-mid` / `fragmented` | The amplicon type that matches your sequencing experiment or protocol. Options are `end-to-end`, `end-to-mid`, or `fragmented`. |
| `--min-coverage` /<br>`-mc`           | Minimum coverage                  | The minimum coverage for the consensus sequence(s). The default is **30**. |
| `--primer-mismatch-rate`/<br>`-pmr`   | Fraction of maximum primer mismatches | The maximum percentage of mismatches allowed during the primer search. Only substitutions are counted (insertions and deletions are not considered).<br>The default is **0.1** (10%). For example, if your primer is 30 nt long, up to 3 mismatches are allowed. |
| `--target` /<br>`--preset`            | Name of the viral target          | A descriptive name for the viral target under analysis (e.g., "Measles", "Sars-CoV-2", "HPV16", or "Influenza_A"). This target will trigger the use of an analysis preset if one is available. If no preset is available, default settings will be used. See [information regarding presets](presets.md).<br>Disable analysis presets using the `--disable-presets` flag. |
| `--disable-presets` /<br>`-dp`        | N/A                               | Disables the use of analysis presets so that default analysis settings are used for all samples or viral targets. The viral target must still be provided. |
| `--match-ref` /<br>`-mr`              | N/A                               | Enables the match-ref process for all samples. See additional information regarding the [match reference process](multi-reference.md#choose-the-best-reference-for-each-sample). |
| `--segmented` /<br>`-seg`             | N/A                               | Indicates that the samples to be analyzed are segmented rather than being based on a single reference genome.<br>This setting applies only to the match-ref process. See details for [segmented viruses](multi-reference.md#choose-the-best-reference-for-each-sample-with-segmented-viruses). |
| `--threads` /<br>`-t`                 | Number of threads                 | The number of local threads available for use. The default is the number of available threads on your system. |
| `--dryrun`                           | N/A                               | Runs the ViroConstrictor workflow without performing any actions. (default: False) |
| `--skip-updates`                     | N/A                               | Skips the check for a new version. |
| `--version` /<br>`-v`                 | N/A                               | Displays the current version of ViroConstrictor and exits. |
| `--help` /<br>`-h`                    | N/A                               | Displays the ViroConstrictor help document and exits. |

## Run a Single-Target Analysis

When running ViroConstrictor as a single-target analysis, you can provide all necessary information via the command line. The provided settings will be applied to all samples in the analysis.

To run an analysis, you need to provide at least the following inputs:

- The input directory containing your FASTQ data via the `--input`/`-i` flag.
- The desired output directory for results and data via the `--output`/`-o` flag.
- The path to a reference FASTA file via the `--reference`/`-ref` flag.
- The sequencing primers in the correct FASTA format (or as a bed-file) via the `--primers`/`-pr` flag.  
    
    !!! Info ""
        If you did not use any primers during your sequencing run, specify `--primers NONE` on the command line.

- A GFF file containing the genomic features for the given reference via the `--features`/`-gff` flag.

    !!! Info ""
        If you do not have access to a GFF file that matches your reference, provide `--features NONE`. In this case, ViroConstrictor will attempt to infer the genomic features during analysis.
        
        However, amino acid translations for the genomic features will **not** be included in the results folder when no GFF file is provided. If amino acid translations are required, supply a matching GFF file for your reference FASTA.

- The sequencing platform via the `--platform` flag; options include "nanopore", "illumina", or "iontorrent".
- The [amplicon type](amplicons.md) in your data via the `--amplicon-type`/`-at` flag.  
    
    !!! Info ""
        If you did not use any primers (and specified `--primers NONE`), this information will be ignored.
        (You still need to provide this flag, even if it is not used.)

- A named viral target via the `--target` or `--preset` flag.

    !!! Attention ""
        The specified viral target will trigger an analysis preset if ViroConstrictor recognizes it as a known target. If no preset is available, default analysis settings will be applied. See [documentation about working with presets](presets.md).  
        You can disable analysis presets with the `--disable-presets` flag.

You may also optionally provide additional information such as the minimum coverage level (`--min-coverage`) and the primer-mismatch rate (`--primer-mismatch-rate`). When these flags are not provided, their default values will be used during analysis.

!!! Example "Start a Single-Target Analysis"
    ```bash
    viroconstrictor \
        --input {path/to/FASTQ-files} \
        --output {path/to/desired/output/folder} \
        --reference {path/to/reference.fasta / NONE} \
        --primers {path/to/primers.fasta / NONE} \
        --features {path/to/features.gff / NONE} \
        --platform {nanopore/illumina/iontorrent} \
        --amplicon-type {end-to-end/end-to-mid} \
        --min-coverage {coverage level} \
        --primer-mismatch-rate {mismatch rate} \
        --target viral-target
    ```

## Run a Multi-Target Analysis

You can run ViroConstrictor in multi-target analysis mode if you wish to set different analysis settings for each sample or target within a single analysis. This is done by providing a samplesheet in Excel, CSV, or TSV format.  
A sample table is shown below, or you can download the example Excel spreadsheet [here](samples.xlsx).

| Sample    | Virus       | Match-ref | Segmented | Primers                              | Reference                                 | Features           | min-coverage | primer-mismatch-rate |
|-----------|-------------|-----------|-----------|--------------------------------------|-------------------------------------------|--------------------|--------------|----------------------|
| Sample_1  | SARS-CoV-2  | FALSE     | FALSE     | /path/to/sars_cov_2_primers.fasta    | /path/to/sars_cov_2_reference.fasta       | /path/to/sars_cov_2.gff | 50           | 0.1                  |
| Sample_2  | Measles     | FALSE     | FALSE     | /path/to/Measles_primers.fasta        | /path/to/measles_reference.fasta          | /path/to/measles.gff | 30           | 0.15                 |
| Sample_3  | Influenza_A | TRUE      | TRUE      | /path/to/Influenza-A_primers.fasta    | /path/to/Influenza-A_reference.fasta       | NONE               | 30           | 0.1                  |
| Sample_4  | HPV         | TRUE      | FALSE     | NONE                                 | /path/to/HPV_reference.fasta              | NONE               | 10           | 0.1                  |

Please note:
- The `Sample` key provided in the samplesheet must correspond exactly to the FASTQ filename in your input directory. If they do not match, ViroConstrictor will alert you to prevent unexpected results.
- The `Virus` column in the samplesheet is used as the viral target name for analysis. This name will trigger an analysis preset if one is recognized. If no preset is available, default analysis settings will be applied. See [documentation about working with presets](presets.md).  
- You can disable analysis presets with the `--disable-presets` flag.

!!! Example "Start a Multi-Target Analysis"
    ```bash
    viroconstrictor \
        --input {path/to/FASTQ-files} \
        --output {path/to/desired/output/folder} \
        --samplesheet {path/to/samplesheet.xlsx} \
        --platform {nanopore/illumina/iontorrent} \
        --amplicon-type {end-to-end/end-to-mid}
    ```

With multi-target analysis, it is also possible to combine run-wide parameters with a samplesheet. For example, if you wish to analyze a batch of samples with different minimum coverage values for each sample, while all other parameters remain the same, you can omit those common details from the samplesheet and specify them on the command line.  
Please note that in the samplesheet, the columns "Sample", "Virus", and "Reference" are mandatory; starting an analysis without these columns will not work.
