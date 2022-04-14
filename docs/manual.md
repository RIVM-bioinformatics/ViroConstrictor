# Starting an analysis

Once ViroConstrictor is [installed](installation.md) you can use it to analyse your NGS FastQ data.  

You can run an analysis in two ways:

1. A [single-target](#run-a-single-target-analysis) analysis with run-wide parameters (the information that you give is applied to every sample in the analysis)
2. A [multi-target](#run-a-multi-target-analysis) analysis where you provide a samplesheet (in excel, .csv or .tsv format) which contains information to run the analysis with different settings for every sample.

Please see `viroconstrictor -h` for all command-line options

---

## Overview of command-line options

Below, you can find a brief summary if all available command line options and their meaning.

| Command | Argument | Explanation |
|---------|----------|-------------|
| `--input` /<br>`-i` | [input directory] | This is the folder containing your input fastq files.<br>Both absolute as relative paths are supported. |
| `--output` /<br>`-o` | [output directory] | This is the desired output folder with all data and results.<br>If the desired folder doesn't exist yet then it will be created. |
| `--reference` /<br>`-ref` | [reference fasta] | Input reference sequence genome in FASTA format |
| `--primers` /<br>`-pr` | [primer fasta] / `NONE` | Used primer sequences in FASTA format.<br>Use `NONE` if your sequencing protocol does not use primers|
| `--features` /<br>`-gff` | [GFF file] / `NONE` | The GFF3 file containing the genomic features (open reading frames etc.) matching the given reference fasta.<br>Use `NONE` if you don't have access to a GFF3 file containing this information.
| `--platform` | `nanopore` / `illumina` / `iontorrent` | The sequencing platform that was used to generate the dataset. Either being 'nanopore', 'illumina' or 'iontorrent'.<br>Default is `nanopore` |
| `--amplicon-type` /<br>`-at` | `end-to-end` / `end-to-mid` | The amplicon type that matches your sequencing experiment/protocol.<br>Either being `end-to-end` or `end-to-mid` |
| `--min-coverage` /<br>`-mc` | Minimum coverage | The minimum coverage for the consensus sequence(s)<br>Default is **30**
| `--primer-mismatch-rate`/<br>`-pmr` | Fraction of maximum primer mismatches | The maximum percentage mismatches that is tolerated during the primer search. Mismatches are counted in substitutions between primer and reference. Insertions and/or deletions are not taken into account.<br>Default is **0.1** (10%)<br>This means that max 10% of the length of a primer may be a mismatch, i.e. if your primer is 30nt then maximum 3 mismatches are allowed |
| `--target` | Name of the viral target | The basic descriptive name of the viral target that is being analysed<br>i.e. "Measles", "Sars-CoV-2", "HPV16", or "Influenza_A"
| `--threads` /<br>`-t` | Amount of threads | Number of local threads that are available to use.<br>Default is the number of available threads in your system |
| `--dryrun` | N/A | Run the ViroConstrictor workflow without actually doing anything.<br>(default: False) |
| `--version` /<br>`-v` | N/A | Shows the current version of ViroConstrictor and exits |
| `--help` /<br>`-h` | N/A | Show the ViroConstrictor help document and exits |


## Run a single-target analysis

When running ViroConstrictor as a single-target analysis its possible to provide all necessary information via the command line. The given information will be applied to all samples in the analysis.

To run an analysis, you need to provide the following inputs/information:  

* The input directory containing your FastQ data via the `--input`/`-i` flag
* The desired output directory for results and data via the `--output`/`-o` flag
* The path to a reference fasta file via the `--reference`/`-ref` flag
* Used sequencing primers [in correct fasta format](#preparing-your-input-primers) or as a bed-file via the `--primers`/`-pr` flag  
    
    !!! Info ""
        If you did not use any primers during your sequencing run then provide `--primers NONE` on the command line

* A GFF file containing the genomic features of the given reference file via the `--features`/`-gff` flag

    !!! Info ""
        If you don't have access to a GFF file with genomic features matching the reference then provide `--features NONE`. ViroConstrictor will then attempt to guess the genomic features of the given reference during analysis.

* The used sequencing platform via `--platform`, can be "nanopore", "illumina", or "iontorrent"
* The [Amplicon-type](#amplicon-types) in your provided data via the `--amplicon-type`/`-at` flag.  
    
    !!! Info ""
        If you did not use any primers, and set `--primers NONE` then this information will be ignored   
        (you still have to provide this command-line flag but it won't be used)

* A named viral target.

Additional information that you can provide, but is not always required is the minimum coverage level (`--min-coverage`) as well as the primer-mismatch rate (`--primer-mismatch-rate`). When these flags are not provided then their default values are used during analysis.


!!! Example "You can start a single-target analysis with a command such as the following:"
    ```bash
    viroconstrictor \
        --input {path/to/FastQ-files} \
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


## Run a multi-target analysis

When running ViroConstrictor for multiple targets in a single analysis. Or if you want to set different analysis settings for every sample then you can run ViroConstrictor in multi-target analysis mode.  
A single analysis with different settings for each sample can be started by providing a samplesheet in Excel, CSV or TSV format.

An example table can be seen below, or download the example excel spreadsheet [here](samples.xlsx)

| Sample | Virus | Match-ref | Primers | Reference | Features | min-coverage | primer-mismatch-rate |
| ------ | ----- | --------- | ------- | --------- | -------- | ------------ | -------------------- |
| Sample_1 | SARS-CoV-2 | FALSE | /path/to/sars_cov_2_primers.fasta | /path/to/sars_cov_2_reference.fasta | /path/to/sars_cov_2.gff | 50 | 0.1 |
| Sample_2 | Measles | TRUE | /path/to/Measles_primers.fasta | /path/to/measles_reference.fasta | /path/to/measles.gff | 30 | 0.15 |
| Sample_3 | Influenza_A | FALSE | /path/to/Influenza-A_primers.fasta | /path/to/Influenza-A_reference.fasta | NONE | 30 | 0.1 |
| Sample_4 | HPV | TRUE | NONE | /path/to/HPV_reference.fasta | NONE | 10 | 0.1 | 

!!! Attention
    Currently the `Match-ref` setting as portrayed in the table above as well as in the example spreadsheet is currently not yet implemented.  
    This is a future functionality that is still being worked on, but we decided to include it already in the samplesheet and command-line arguments as we don't want to significantly change these formats going forward (unless actually deemed necessary).  
    Setting the `Match-ref` to `TRUE` as in the example above will therefore not do anything in the analysis as per the current version.

Please keep in mind that the Sample key as given in de samplesheet must correspond with the FastQ filename in your input directory. If this isn't the case then ViroConstrictor will let you know without running to ensure you won't get unexplainable data.

!!! Example "You can start a multi-target analysis with a command such as the following"
    ```bash
    viroconstrictor \
        --input {path/to/FastQ-files} \
        --output {path/to/desired/output/folder} \
        --samplesheet {path/to/samplesheet.xlsx} \
        --platform {nanopore/illumina/iontorrent} \
        --amplicon-type {end-to-end/end-to-mid}
    ```

With a multi-target analysis it's also possible to combine run-wide parameters with a samplesheet.  
For example, you wish to analyse a batch of samples and set different minimum coverage values for every sample, but all other parameters in your analysis are the same.  
You can then leave out certain information in your samplesheet and supplement this information through the command line.  
Please note that in the given samplesheet, the columns "Samples", "Virus", and "Reference" are considered to be mandatory. starting an analysis without this information in the samplesheet will not work.
