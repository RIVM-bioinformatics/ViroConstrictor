# User Guide

Once ViroConstrictor is [installed](installation.md) you can use it to analyse your NGS FastQ data.  

To run an analysis, you need the following inputs:  

* Input FASTQ data
* Reference fasta
* Used sequencing primers [in correct format](#preparing-your-input-primers)
* Sequencing platform
* [Amplicon-type](#amplicon-types)

## Preparing your input primers

In order to get optimal results, please make sure the fasta headers in your fasta file with primers are formatted properly. Please make sure the fasta headers are formatted according the following format:

`>{primer-name}_{primer-number}_{orientation}`

It's important that the primers which together form a single amplicon have the same primer-name and number.  

Orientation keywords for forward primers are: *"LEFT"*/*"PLUS"*/*"POSITIVE"*/*"FORWARD"*  
Orientation keywords for reverse primers are: *"RIGHT"*/*"MINUS"*/*"NEGATIVE"*/*"REVERSE"*

!!! example "Example of formatted primer names from the ArticV3 SARS-CoV-2 sequencing protocol"
    ```Markdown
    >nCoV-2019_1_LEFT  
    ACCAACCAACTTTCGATCTCTTGT  
    >nCoV-2019_1_RIGHT  
    CATCTTTAAGATGTTGACGTGCCTC  
    >nCoV-2019_2_LEFT  
    CTGTTTTACAGGTTCGCGACGT  
    >nCoV-2019_2_RIGHT  
    TAAGGATCAGTGCCAAGCTCGT
    ```

If your protocol has alternative primers then make sure the fasta header contains the "alt" keyword in the following format:

`>{primer-name}_{primer-number}_alt_{orientation}`  

Please make sure the "alt" keyword is in the middle and not at the end of the fasta header.

!!! example "Example of formatted primer names from the ArticV3 SARS-CoV-2 sequencing protocol with alternative primers included"
    ```Markdown
    >nCoV-2019_13_LEFT  
    TCGCACAAATGTCTACTTAGCTGT  
    >nCoV-2019_13_RIGHT  
    ACCACAGCAGTTAAAACACCCT  
    >nCoV-2019_14_LEFT  
    CATCCAGATTCTGCCACTCTTGT  
    >nCoV-2019_14_alt_LEFT  
    TGGCAATCTTCATCCAGATTCTGC  
    >nCoV-2019_14_RIGHT  
    AGTTTCCACACAGACAGGCATT  
    >nCoV-2019_14_alt_RIGHT  
    TGCGTGTTTCTTCTGCATGTGC  
    >nCoV-2019_15_LEFT  
    ACAGTGCTTAAAAAGTGTAAAAGTGCC  
    >nCoV-2019_15_alt_LEFT  
    AGTGCTTAAAAAGTGTAAAAGTGCCT  
    >nCoV-2019_15_RIGHT  
    AACAGAAACTGTAGCTGGCACT  
    >nCoV-2019_15_alt_RIGHT  
    ACTGTAGCTGGCACTTTGAGAGA
    ```  

## Amplicon types

The "amplicon type" is basically a description of a PCR amplicon, and which part of said amplicon is actually sequenced.  
This information is important because given primers will be removed from the reads during the analysis.

In the ViroConstrictor pipeline the primer-removal is done by our tool [AmpliGone](https://rivm-bioinformatics.github.io/AmpliGone/).

The information regarding the "amplicon type" is required because it tells te pipeline, and therefore also AmpliGone, from which side(s) of a read the primer-sequence has to be removed.

Currently there are two options: "End-to-End" or "End-to-Mid".  

In short: "End-to-End" means that the sequenced read covers the full length of an amplicon. Meaning that the primer-sequence is present at both ends of a read.  
"End-to-Mid" means that the sequences read only *partially* covers the length of the amplicon. Meaning that the primer-sequence is present at only one end of a read.

Also see the [AmpliGone documentation regarding amplicon-types](https://rivm-bioinformatics.github.io/AmpliGone/latest/amplicon-types/) where we go into more detail regarding these types and what these terms mean.  
Please check your sequencing and laboratory setup for the amplicon type to match your analysis to ensure the best results.

## Running an analysis

Please also see the command-line help document for a quick explanation of every possible argument: `viroconstrictor -h`

!!! example "You can start an analysis with a command such as the following:"
    ```bash
    viroconstrictor \
        --input {path/to/FastQ-files} \
        --output {path/to/desired/output/folder} \
        --primers {path/to/primers.fasta / NONE} \
        --platform {nanopore/illumina/iontorrent} \
        --amplicon-type {end-to-end/end-to-mid} \
        --threads {threads}
    ```

Here, `threads` refers to the amount of threads that ViroConstrictor may use on your LOCAL machine. If you're using ViroConstrictor on a HPC/cluster then these threads will only be used during the pre-processing steps.  
If you're using ViroConstrictor on your local machine then the given amount of threads will act as a 'ceiling' of usable threads during analysis.


If your sequencing protocol does not use primers then set "NONE" for the primers flag with `--primers NONE`.

Below, you can find a brief summary if all available command line options and their meaning.

| Command | Argument | Explanation |
|---------|----------|-------------|
| `--input` /<br>`-i` | [input directory] | This is the folder containing your input fastq files.<br>Both absolute as relative paths are supported. |
| `--output` /<br>`-o` | [output directory] | This is the desired output folder with all data and results.<br>If the desired folder doesn't exist yet then it will be created. |
| `--reference` /<br>`-ref` | [reference fasta] | Input reference sequence genome in FASTA format |
| `--primers` /<br>`-pr` | [primer fasta] / `NONE` | Used primer sequences in FASTA format.<br>Use `NONE` if your sequencing protocol does not use primers|
| `--platform` | `nanopore` / `illumina` / `iontorrent` | The sequencing platform that was used to generate the dataset. Either being 'nanopore', 'illumina' or 'iontorrent'.<br>Default is `nanopore` |
| `--amplicon-type` /<br>`-at` | `end-to-end` / `end-to-mid` | The amplicon type that matches your sequencing experiment/protocol.<br>Either being `end-to-end` or `end-to-mid` |
| `--primer-mismatch-rate`/<br>`-pmr` | number of maximum mismatches | The maximum amount of mismatches that are tolerated during the primer search. Mismatches are counted in substitutions between primer and reference. Insertions and/or deletions are not taken into account.<br>Default maximum amount of mismatches is **3** |
| `--features` /<br>`-gff` | [GFF file] / `NONE` | A GFF3 file that contains the Open Reading Frame information of the given reference.<br>Use `NONE` if you don't have access to a GFF3 file with this information |
| `--threads` /<br>`-t` | Amount of threads | Number of local threads that are available to use.<br>Default is the number of available threads in your system |
| `--dryrun` | N/A | Run the ViroConstrictor workflow without actually doing anything.<br>(default: False) |
| `--version` /<br>`-v` | N/A | Shows the current version of ViroConstrictor and exits |
| `--help` /<br>`-h` | N/A | Show the ViroConstrictor help document and exits |