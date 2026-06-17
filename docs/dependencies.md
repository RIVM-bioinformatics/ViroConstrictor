# Key Dependencies

ViroConstrictor relies on several important libraries and tools to provide robust viral genome analysis.

## System Requirements

ViroConstrictor requires Python 3.10 or higher and runs on Linux systems. macOS may work but is not extensively tested. Windows is not supported (use WSL or containers instead).

## Core Dependencies

**[Snakemake](https://snakemake.readthedocs.io/)** - Workflow management and pipeline orchestration with parallel processing support

**[Biopython](https://biopython.org/)** - Biological sequence analysis, FASTA/GenBank parsing, and format conversions

**[BioValid](https://github.com/RIVM-bioinformatics/biovalid)** - Biological data validation for input files and sequences

**[Pandas](https://pandas.pydata.org/)** - Data manipulation, samplesheet processing, and result aggregation

**[PyYAML](https://pyyaml.org/wiki/PyYAMLDocumentation)** - Configuration file handling and workflow parameter management

**[OpenPyXL](https://openpyxl.readthedocs.io/)** - Excel samplesheet reading and writing support

**[Rich](https://rich.readthedocs.io/)** - Enhanced terminal output with colored formatting and progress bars

**[DRMAA](https://drmaa-python.readthedocs.io/en/latest/)** - Grid computing integration for HPC cluster job submission (SLURM, LSF)

## Workflow Tools

**[Minimap2](https://github.com/lh3/minimap2)** - Fast sequence alignment for mapping reads to reference genomes

**[Samtools](https://www.htslib.org/)** - SAM/BAM file manipulation and processing

**[Pysam](https://pysam.readthedocs.io/)** - Python interface for SAM/BAM file reading and writing

**[FastP](https://github.com/OpenGene/fastp)** - Quality control and adapter trimming for sequencing data

**[AmpliGone](https://rivm-bioinformatics.github.io/AmpliGone/)** - Primer removal from amplicon-based sequencing data

**[TrueConsense](https://rivm-bioinformatics.github.io/TrueConsense/latest/)** - Consensus sequence generation from aligned reads

**[AminoExtract](https://github.com/RIVM-bioinformatics/AminoExtract)** - Amino acid sequence extraction and translation from genomic features

**[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)** - Quality assessment of sequencing data

**[MultiQC](https://multiqc.info/)** - Aggregated quality control reporting

**[BEDtools](https://bedtools.readthedocs.io/)** - Genomic interval operations and analysis