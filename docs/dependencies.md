# Key Dependencies

ViroConstrictor relies on several important libraries and tools to provide robust viral genome analysis.

## System Requirements

ViroConstrictor requires Python 3.10 or higher and runs on Linux systems. macOS may work but is not extensively tested. Windows is not supported (use WSL or containers instead).

## Core Dependencies

**[Snakemake](https://snakemake.readthedocs.io/)** [^1] - Workflow management and pipeline orchestration with parallel processing support

**[Biopython](https://biopython.org/)** [^2] - Biological sequence analysis, FASTA/GenBank parsing, and format conversions

**[BioValid](https://github.com/RIVM-bioinformatics/biovalid)** [^3] - Biological data validation for input files and sequences

**[Pandas](https://pandas.pydata.org/)** [^4] - Data manipulation, samplesheet processing, and result aggregation

**[PyYAML](https://pyyaml.org/wiki/PyYAMLDocumentation)** [^5] - Configuration file handling and workflow parameter management

**[OpenPyXL](https://openpyxl.readthedocs.io/)** [^6] - Excel samplesheet reading and writing support

**[Rich](https://rich.readthedocs.io/)** [^7] - Enhanced terminal output with colored formatting and progress bars

**[DRMAA](https://drmaa-python.readthedocs.io/en/latest/)** [^8] - Grid computing integration for HPC cluster job submission (SLURM, LSF)

## Workflow Tools

**[Minimap2](https://github.com/lh3/minimap2)** [^9] - Fast sequence alignment for mapping reads to reference genomes

**[Samtools](https://www.htslib.org/)** [^10] - SAM/BAM file manipulation and processing

**[Pysam](https://pysam.readthedocs.io/)** [^11] - Python interface for SAM/BAM file reading and writing

**[FastP](https://github.com/OpenGene/fastp)** [^12] - Quality control and adapter trimming for sequencing data

**[AmpliGone](https://rivm-bioinformatics.github.io/AmpliGone/)** [^13] - Primer removal from amplicon-based sequencing data

**[TrueConsense](https://rivm-bioinformatics.github.io/TrueConsense/latest/)** [^14] - Consensus sequence generation from aligned reads

**[AminoExtract](https://github.com/RIVM-bioinformatics/AminoExtract)** [^15] - Amino acid sequence extraction and translation from genomic features

**[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)** [^16] - Quality assessment of sequencing data

**[MultiQC](https://multiqc.info/)** [^17] - Aggregated quality control reporting

**[BEDtools](https://bedtools.readthedocs.io/)** [^18] - Genomic interval operations and analysis






[^1]:  Mölder F, Jablonski KP, Letcher B et al. Sustainable data analysis with Snakemake [version 3; peer review: 2 approved]. F1000Research 2025, 10:33 (https://doi.org/10.12688/f1000research.29032.3) 
[^2]: Cock PJA, et al. Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics. 2009;25(11):1422–1423 (https://doi.org/10.1093/bioinformatics/btp163)
[^3]: RIVM-bioinformatics. BioValid Documentation (https://github.com/RIVM-bioinformatics/biovalid)
[^4]: McKinney W. Data Structures for Statistical Computing in Python. Proceedings of the 9th Python in Science Conference. 2010;445:51–56 (https://doi.org/10.25080/Majora-92bf1922-00a)
[^5]: Kirill Simonov. PyYAML Documentation (https://pyyaml.org/wiki/PyYAMLDocumentation)
[^6]: openpyxl Documentation (https://openpyxl.readthedocs.io/)
[^7]: Will McGugan. Rich Documentation (https://rich.readthedocs.io/)
[^8]: DRMAA Python Documentation (https://drmaa-python.readthedocs.io/en/latest/)
[^9]: Li H. Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics. 2018;34(18):3094–3100 (https://doi.org/10.1093/bioinformatics/bty191)
[^10]: Danecek P, et al. Twelve years of SAMtools and BCFtools. GigaScience. 2021;10(2):giab008 (https://doi.org/10.1093/gigascience/giab008)
[^11]: Pysam Documentation (https://pysam.readthedocs.io/)
[^12]: Chen S, Zhou Y, Chen Y, Gu J. fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics. 2018;34(17):i884–i890 (https://doi.org/10.1093/bioinformatics/bty560)
[^13]: AmpliGone Documentation (https://rivm-bioinformatics.github.io/AmpliGone/)
[^14]: TrueConsense Documentation (https://rivm-bioinformatics.github.io/TrueConsense/latest/)
[^15]: AminoExtract Documentation (https://github.com/RIVM-bioinformatics/AminoExtract)
[^16]: Andrews S. FastQC: a quality control tool for high throughput sequence data. 2010 (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
[^17]: Ewels P, Magnusson M, Lundin S, Käller M. MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics. 2016;32(19):3047–3048 (https://doi.org/10.1093/bioinformatics/btw354)
[^18]: Quinlan AR, Hall IM. BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics. 2010;26(6):841–842 (https://doi.org/10.1093/bioinformatics/btq033)