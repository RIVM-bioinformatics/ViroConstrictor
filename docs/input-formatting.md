# Formatting your inputs

For optimal results and accurate data analysis, it is important that the various inputs you provide are formatted correctly. Input formatting is particularly crucial for the primer file, the reference FASTA, and the GFF file.

Please note that ViroConstrictor does not support circular reference data. Analysing a circular virus, such as HPV, is possible as long as the analysis uses a linear input reference with matching primers and GFF.

---

## Formatting your input primers

ViroConstrictor accepts two primer input formats: [BED](https://en.wikipedia.org/wiki/BED_(file_format)) and [FASTA](https://en.wikipedia.org/wiki/FASTA_format). The formatting of primer names directly impacts your analysis results.

### Primers in FASTA format

For optimal results, ensure that the FASTA headers in your primer file are formatted correctly.

The FASTA headers for your primers should follow this format:  
`>{primer-name}_{primer-number}_{orientation}`

It is important that primers forming a single amplicon share the same primer name and number.

Orientation keywords for forward primers include: *"LEFT"*, *"PLUS"*, *"POSITIVE"*, and *"FORWARD"*.  
Orientation keywords for reverse primers include: *"RIGHT"*, *"MINUS"*, *"NEGATIVE"*, and *"REVERSE"*.

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

If your protocol includes alternative primers, ensure the FASTA header contains the "alt" keyword in the following format:

`>{primer-name}_{primer-number}_alt_{orientation}`  

The "alt" keyword must appear in the middle and not at the end of the header.

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

### Primers in BED format

As with FASTA, it is important that the primer name in the BED file is formatted properly to determine which primers form a single amplicon. The 'strand' column is used to determine the orientation of the primer.

Please ensure the following:

* The first column (the reference ID column) must contain the identifier that matches the identifier in your reference FASTA.

* The fourth column, which contains the primer name, should follow this format:  
  `{primer-name}_{primer-number}_{orientation}`  
    
    !!! info ""
        If your protocol includes alternative primers, ensure the primer name includes the "alt" keyword as follows:  
        `{primer-name}_{primer-number}_alt_{orientation}`

* The sixth column (the strand column) indicates the primer orientation:  
  `+` for forward and `-` for reverse, relative to the reference.

The table below illustrates what your BED file with primer information should look like (note that a real BED file does not include headers):

| Reference ID | start coordinate | stop coordinate | primer name        | score | strand |
|--------------|------------------|-----------------|--------------------|-------|--------|
| MN908947.3   | 25               | 50              | ncov-2019_1_LEFT   | .     | +      |
| MN908947.3   | 324              | 344             | ncov-2019_2_LEFT   | .     | +      |
| MN908947.3   | 408              | 431             | ncov-2019_1_RIGHT  | .     | -      |
| MN908947.3   | 644              | 666             | ncov-2019_3_LEFT   | .     | +      |
| MN908947.3   | 705              | 727             | ncov-2019_2_RIGHT  | .     | -      |

## Formatting your input reference FASTA

The formatting of your reference FASTA is straightforward, with no strict formatting requirements. However, ensure that your reference sequence has a proper identifier in the FASTA header, as this identifier is used during analysis.

Nucleotide ambiguity codes in your reference FASTA are supported but generally discouraged, as they can negatively affect the primer removal process. You will receive a warning during pre-processing if ambiguous nucleotides are detected in your reference sequence.

If your reference FASTA contains multiple sequences, analysis will be performed for each individual sequence, split based on the reference identifier. This feature can be useful for analysing segmented viruses such as Influenza.

## Formatting your input GFF

Formatting the input GFF is usually optional and minimal.

If you want proper extraction and translation of amino acids from the generated consensus sequence, ensure that the input GFF includes either a "Name" or an "ID" attribute in the attributes column.

For example: `Name="Nucleocapsid"` or `ID="N"` in the GFF attributes column.

This information is used to group the extracted amino acid sequences from samples during analysis.
