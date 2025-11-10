# Quick Run Example

## Install

```bash
# Install ViroConstrictor with Mamba
mamba create --name viroconstrictor -c conda-forge -c bioconda viroconstrictor
mamba activate viroconstrictor
```

---

## Run Analysis

```bash
viroconstrictor \
    --input /path/to/FASTQ-files \
    --output /path/to/output \
    --reference /path/to/reference.fasta \
    --primers /path/to/primers.fasta \
    --features /path/to/features.gff \
    --platform nanopore \
    --amplicon-type end-to-end \
    --target viral-target
```

---

For comprehensive installation options and detailed usage, see the [full installation guide](installation.md) and [manual](manual.md).