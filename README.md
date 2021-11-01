# ViroConstrictor

ViroConstrictor is a pipeline designed to process raw FastQ data from viral amplicon-based sequencing and generate  biologically correct consensus sequences from your data based on a given reference genome.

ViroConstrictor performs high speed data quality control, data cleanup and high accuracy removal of primer-sequences from NGS reads. As well as alignment of reads and generation of a consensus sequence using the TrueConsense consensus-caller which accounts for sequencing errors and alignment artefacts.

ViroConstrictor is able to run both on a standalone (linux) computer, as well as High-Performance Computing (HPC) infrastructures. 

ViroConstrictor is compatible with Nanopore, Illumina, and IONTorrent data (fastq).

Please see [the documentation](https://rivm-bioinformatics.github.io/ViroConstrictor/) for more information.

ViroConstrictor is available under the [AGPLv3 licence](https://www.gnu.org/licenses/agpl-3.0.en.html)