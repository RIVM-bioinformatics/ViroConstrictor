# Key Dependencies

ViroConstrictor relies on several important libraries and tools to provide robust viral genome analysis.

## System Requirements

ViroConstrictor requires Python 3.10 or higher and runs on Linux systems. macOS may work but is not extensively tested. Windows is not supported (use WSL or containers instead).

## Core Dependencies

**Snakemake**
- Workflow management and pipeline orchestration with parallel processing support

**Biopython** 
- Biological sequence analysis, FASTA/GenBank parsing, and format conversions

**BioValid**
- Biological data validation for input files and sequences

**Pandas**
- Data manipulation, samplesheet processing, and result aggregation

**PyYAML**
- Configuration file handling and workflow parameter management

**OpenPyXL**
- Excel samplesheet reading and writing support

**Rich**
- Enhanced terminal output with colored formatting and progress bars

**DRMAA**
- Grid computing integration for HPC cluster job submission (SLURM, LSF)

## Workflow Tools

**Minimap2**
- Fast sequence alignment for mapping reads to reference genomes

**Samtools**
- SAM/BAM file manipulation and processing

**Pysam**
- Python interface for SAM/BAM file reading and writing

**FastP**
- Quality control and adapter trimming for sequencing data

**AmpliGone**
- Primer removal from amplicon-based sequencing data

**TrueConsense**
- Consensus sequence generation from aligned reads

**AminoExtract**
- Amino acid sequence extraction and translation from genomic features

**FastQC**
- Quality assessment of sequencing data

**MultiQC**
- Aggregated quality control reporting

**BEDtools**
- Genomic interval operations and analysis