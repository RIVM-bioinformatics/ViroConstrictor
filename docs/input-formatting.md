# Formatting Your Input Files

For optimal results and accurate analysis, it is important that your input files are formatted correctly. This guide covers the formatting requirements for primers, reference sequences, and annotation files.

!!! warning "Circular Reference Limitation"
    ViroConstrictor does not support circular reference data. To analyze circular viruses (such as HPV), use a linear input reference with matching primers and annotations.

---

## Input File Types Overview

ViroConstrictor accepts several input file formats:

- **Primers**: FASTA or BED format
- **Reference sequences**: FASTA or GenBank format  
- **Genomic features**: GFF3 format (optional when using GenBank)

## Reference Sequences

### FASTA Format

The formatting of your reference FASTA is straightforward, with no strict formatting requirements. However, ensure that your reference sequence has a proper identifier in the FASTA header, as this identifier is used during analysis.

<!-- TODO: Add example of good vs bad FASTA headers -->

Nucleotide ambiguity codes in your reference FASTA are supported but generally discouraged, as they can negatively affect the primer removal process. You will receive a warning during pre-processing if ambiguous nucleotides are detected in your reference sequence.

If your reference FASTA contains multiple sequences, analysis will be performed for each individual sequence, split based on the reference identifier. This feature can be useful for analysing segmented viruses such as Influenza.

### GenBank Format

ViroConstrictor also accepts GenBank files (`.gb` or `.gbk`) as reference input. GenBank files contain both the reference sequence and genomic feature annotations in a single file, which offers several advantages:

- **Simplified workflow**: No need to provide separate GFF files
- **Integrated annotations**: Features like genes, CDS, and other annotations are included
- **Automatic feature extraction**: ViroConstrictor automatically extracts genomic features for amino acid translation
- **No viral target required**: When using GenBank format, you don't need to specify a viral target

When using a GenBank file, simply provide it via the `--reference` parameter, and set `--features NONE` since the annotations are already included in the GenBank file.

<!-- TODO: Consider adding guidance on where to find quality GenBank files (NCBI, etc.) -->

---

## Primer Files

ViroConstrictor accepts two primer input formats: [BED](https://en.wikipedia.org/wiki/BED_(file_format)) and [FASTA](https://en.wikipedia.org/wiki/FASTA_format). The formatting of primer names directly impacts your analysis results.

<!-- TODO: Explain what happens if primer names are formatted incorrectly -->

### FASTA Format

For optimal results, ensure that the FASTA headers in your primer file are formatted correctly.

**Required format**: `>{primer-name}_{primer-number}_{orientation}`

**Key requirements**:
- Primers forming a single amplicon must share the same primer name and number
- Use consistent orientation keywords

**Accepted orientation keywords**:
- **Forward primers**: "LEFT", "PLUS", "POSITIVE", "FORWARD"
- **Reverse primers**: "RIGHT", "MINUS", "NEGATIVE", "REVERSE"

<!-- TODO: Clarify whether orientation keywords are case-sensitive -->

!!! example "Example of formatted primer names from the ArticV3 SARS-CoV-2 sequencing protocol"
    ```markdown
    >nCoV-2019_1_LEFT  
    ACCAACCAACTTTCGATCTCTTGT  
    >nCoV-2019_1_RIGHT  
    CATCTTTAAGATGTTGACGTGCCTC  
    >nCoV-2019_2_LEFT  
    CTGTTTTACAGGTTCGCGACGT  
    >nCoV-2019_2_RIGHT  
    TAAGGATCAGTGCCAAGCTCGT
    ```

**Alternative primers format**: `>{primer-name}_{primer-number}_alt_{orientation}`

!!! important "Alternative Primer Placement"
    The "alt" keyword must appear in the middle of the header, not at the end.

!!! example "Example of formatted primer names from the ArticV3 SARS-CoV-2 sequencing protocol with alternative primers included"
    ```markdown
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

### BED Format

BED format provides coordinate-based primer information. The primer name formatting follows the same rules as FASTA format.

**Key requirements**:
- **Column 1** (reference ID): Must match the identifier in your reference FASTA
- **Column 4** (primer name): Use format `{primer-name}_{primer-number}_{orientation}` (or `{primer-name}_{primer-number}_alt_{orientation}` for alternatives)
- **Column 6** (strand): `+` for forward primers, `-` for reverse primers

<!-- TODO: Add guidance on coordinate accuracy and 0-based vs 1-based indexing -->

The table below illustrates what your BED file with primer information should look like (note that a real BED file does not include headers):

| Reference ID | start coordinate | stop coordinate | primer name        | score | strand |
|--------------|------------------|-----------------|--------------------|-------|--------|
| MN908947.3   | 25               | 50              | ncov-2019_1_LEFT   | .     | +      |
| MN908947.3   | 324              | 344             | ncov-2019_2_LEFT   | .     | +      |
| MN908947.3   | 408              | 431             | ncov-2019_1_RIGHT  | .     | -      |
| MN908947.3   | 644              | 666             | ncov-2019_3_LEFT   | .     | +      |
| MN908947.3   | 705              | 727             | ncov-2019_2_RIGHT  | .     | -      |

---

## Genomic Features (GFF3 Files)

GFF3 files provide genomic feature annotations and are optional when using GenBank references.

**When to use GFF3 files**:
- Using FASTA reference files
- Need amino acid translations from consensus sequences
- Want custom feature annotations

**Formatting requirements**:
- Include either "Name" or "ID" attribute in the attributes column
- Example: `Name="Nucleocapsid"` or `ID="N"`
- This information groups extracted amino acid sequences during analysis

<!-- TODO: Add examples of common GFF3 formatting issues and how to fix them -->

**When GFF3 files are not needed**:
- Using GenBank reference files (annotations included)
- Only interested in nucleotide consensus sequences
- Using analysis presets that don't require amino acid translations
