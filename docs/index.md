---
hide:
  - navigation
  - toc
---

[![CodeFactor](https://www.codefactor.io/repository/github/rivm-bioinformatics/viroconstrictor/badge)](https://www.codefactor.io/repository/github/rivm-bioinformatics/viroconstrictor)
![Snakemake](https://img.shields.io/badge/snakemake-6.4.1-brightgreen.svg?style=flat-square)

![GitHub release (latest by date including pre-releases)](https://img.shields.io/github/v/release/RIVM-bioinformatics/ViroConstrictor?include_prereleases)
![GitHub](https://img.shields.io/github/license/RIVM-bioinformatics/ViroConstrictor)

# ViroConstrictor

ViroConstrictor is a pipeline designed to process raw FastQ data from amplicon-based sequencing experiments and generate a biologically correct consensus sequence of the sequenced genome.

ViroConstrictor performs high speed data quality control, data cleanup and high accuracy removal of primer-sequences from NGS reads. As well as alignment of reads and generation of a consensus sequence using the TrueConsense consensus-caller which accounts for sequencing errors and alignment artefacts.

ViroConstrictor is able to run both on a standalone (linux) computer, as well as High-Performance Computing (HPC) infrastructures.

ViroConstrictor is compatible with Nanopore, Illumina, and IONTorrent data (fastq).

Please see the [installation page](installation.md) for download and installation instructions, and please see the [user guide page](user-guide.md) page for basic usage instructions.

ViroConstrictor is available under the [AGPLv3 licence](https://www.gnu.org/licenses/agpl-3.0.en.html) 

---
## Authors

* Florian Zwagemaker
* Dennis Schmitz
* Karim Hajji
* Annelies Kroneman