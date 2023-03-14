# ViroConstrictor

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/viroconstrictor/README.html)

![GitHub release (latest by date including pre-releases)](https://img.shields.io/github/v/release/RIVM-bioinformatics/ViroConstrictor?include_prereleases)
![Conda](https://img.shields.io/conda/v/bioconda/viroconstrictor)  
![GitHub](https://img.shields.io/github/license/RIVM-bioinformatics/ViroConstrictor)
[![CodeFactor](https://www.codefactor.io/repository/github/rivm-bioinformatics/viroconstrictor/badge)](https://www.codefactor.io/repository/github/rivm-bioinformatics/viroconstrictor)  
![GitHub deployments](https://img.shields.io/github/deployments/RIVM-bioinformatics/ViroConstrictor/github-pages?label=Documentation%20deployment)  

ViroConstrictor is a pipeline designed to process raw FastQ data from viral amplicon-based sequencing and generate  biologically correct consensus sequences from your data based on a given reference genome.

ViroConstrictor performs high speed data quality control, data cleanup and high accuracy removal of primer-sequences from NGS reads. As well as alignment of reads and generation of a consensus sequence using the TrueConsense consensus-caller which accounts for sequencing errors and alignment artefacts.

ViroConstrictor is able to run both on a standalone (linux) computer, as well as High-Performance Computing (HPC) infrastructures. 

ViroConstrictor is compatible with Nanopore, Illumina, and IONTorrent data (fastq).

Please see [the documentation](https://rivm-bioinformatics.github.io/ViroConstrictor/) for more information.

ViroConstrictor is available under the [AGPLv3 licence](https://www.gnu.org/licenses/agpl-3.0.en.html)