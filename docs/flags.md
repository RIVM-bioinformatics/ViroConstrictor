
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
