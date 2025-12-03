---
hide:
  - navigation
  - toc
---
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/viroconstrictor/README.html)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7688035.svg)](https://doi.org/10.5281/zenodo.7688035)

[![GitHub release (latest by date including pre-releases)](https://img.shields.io/github/v/release/RIVM-bioinformatics/ViroConstrictor?include_prereleases)](https://github.com/RIVM-bioinformatics/ViroConstrictor/releases/latest)
[![Conda](https://img.shields.io/conda/v/bioconda/viroconstrictor)](https://anaconda.org/bioconda/viroconstrictor)  
[![GitHub](https://img.shields.io/github/license/RIVM-bioinformatics/ViroConstrictor)](https://github.com/RIVM-bioinformatics/ViroConstrictor/blob/main/LICENSE)
[![CodeFactor](https://www.codefactor.io/repository/github/rivm-bioinformatics/viroconstrictor/badge)](https://www.codefactor.io/repository/github/rivm-bioinformatics/viroconstrictor)  


# ViroConstrictor

ViroConstrictor is a pipeline designed to process raw FASTQ data from amplicon-based sequencing experiments and generate a biologically accurate consensus sequence of the sequenced genome.

ViroConstrictor performs high-speed data quality control, data cleanup, and highly accurate removal of primer sequences from NGS reads, as well as aligning reads and generating a consensus sequence using the [TrueConsense](https://rivm-bioinformatics.github.io/TrueConsense/0.5.1/) consensus caller, which accounts for sequencing errors and alignment artefacts.

ViroConstrictor can run on both standalone (Linux) computers and high-performance computing (HPC) infrastructures.

ViroConstrictor is compatible with Nanopore, Illumina, and IonTorrent data (FASTQ).

Please see the [getting started](installation.md) page for download and installation instructions, and refer to the [user guide page](manual.md) for basic usage instructions.

ViroConstrictor is available under the [AGPLv3 licence](https://www.gnu.org/licenses/agpl-3.0.en.html).

---
## Authors

* Florian Zwagemaker
* Dennis Schmitz
* Karim Hajji
* Gino Raaijmakers
* Iveta Krskova
* Annelies Kroneman
* Margo Raijmakers


If you use ViroConstrictor in your work, please cite:
> Zwagemaker, F., Hajji, K., Raaijmakers, G., Schmitz, D., Raijmakers, M., Kršková, I., Kroneman, A., & The RIVM-IDS Bioinformatics team. ViroConstrictor [Computer software]. https://doi.org/10.5281/zenodo.7688035 