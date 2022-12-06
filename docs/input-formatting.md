# Formatting your inputs

For optimal results and the best analysis of your data, it is important that the various inputs you may provide are formatted correctly.  
Formatting of the input is particularly of importance for the input primers file, the reference fasta, and the GFF file.

When formatting the inputs please note that ViroConstrictor does not support circular reference information (and related information). Analysing a circular virus, such as HPV, is possible as long as the analysis itself happens with a 'linear' input reference and matching primers/GFF.

---

## Formatting your input primers

ViroConstrictor has two options for primer input: [BED](https://en.wikipedia.org/wiki/BED_(file_format)) and [fasta](https://en.wikipedia.org/wiki/FASTA_format).  
The formatting of the input primer names may directly impact your analysis results.

### Primers in Fasta format

In order to get optimal results, please make sure the fasta headers in your fasta file with primers are formatted properly.  

Please make sure the fasta headers for your primers are formatted according the following format:  
`>{primer-name}_{primer-number}_{orientation}`

It's important that the primers which together form a single amplicon have the same primer-name and number.  

Orientation keywords for forward primers are: *"LEFT"*/*"PLUS"*/*"POSITIVE"*/*"FORWARD"*  
Orientation keywords for reverse primers are: *"RIGHT"*/*"MINUS"*/*"NEGATIVE"*/*"REVERSE"*

!!! example "Example of formatted primer names from the ArticV3 SARS-CoV-2 sequencing protocol"
    ```Markdown
    >nCoV-2019_1_LEFT  
    ACCAACCAACTTTCGATCTCTTGT  
    >nCoV-2019_1_RIGHT  
    CATCTTTAAGATGTTGACGTGCCTC  
    >nCoV-2019_2_LEFT  
    CTGTTTTACAGGTTCGCGACGT  
    >nCoV-2019_2_RIGHT  
    TAAGGATCAGTGCCAAGCTCGT
    ```

If your protocol has alternative primers then make sure the fasta header contains the "alt" keyword in the following format:

`>{primer-name}_{primer-number}_alt_{orientation}`  

Please make sure the "alt" keyword is in the middle and not at the end of the fasta header.

!!! example "Example of formatted primer names from the ArticV3 SARS-CoV-2 sequencing protocol with alternative primers included"
    ```Markdown
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

As with the primers in fasta format, it is important the primer name is formatted properly in the given BED file. This is mainly important for determining which primers together form a single amplicon.  
When using a BED file, the 'strand' column is leading in determining the orientation of the primer.

Please make sure of the following:

* The content of first column, described as the reference ID column underneath, must contain the identifier matching the identifier of your reference fasta.

* The fourth column, being the primer name, follows the following format:  
`{primer-name}_{primer-number}_{orientation}`  
    
    !!! info ""
        if your protocol has alternative primers then make sure the primer name contains the "alt" keyword as follows:  
        `{primer-name}_{primer-number}_alt_{orientation}`

* The sixth colum, being the "strand" column, indicates the orientation of the primer.  
`+` being forward and `-` being reverse, relative to the reference.

The table below gives you an impression what your BED file with primer information should look like, please note however that a real BED file does not have headers. 

| Reference ID 	| start coordinate  	| stop coordinate  	| primer name  	        | score | strand    |
|------------	|-----	                |-----	            |-------------------	|----	|-------    |
| MN908947.3 	| 25  	                | 50  	            | ncov-2019_1_LEFT  	| . 	| + 	    |
| MN908947.3 	| 324 	                | 344 	            | ncov-2019_2_LEFT  	| . 	| + 	    |
| MN908947.3 	| 408 	                | 431 	            | ncov-2019_1_RIGHT 	| . 	| - 	    |
| MN908947.3 	| 644 	                | 666 	            | ncov-2019_3_LEFT  	| . 	| + 	    |
| MN908947.3 	| 705 	                | 727 	            | ncov-2019_2_RIGHT 	| . 	| - 	    |


## Formatting your input reference fasta

The formatting of your reference fasta is pretty straightforward, and there are no breaking requirements.  

Make sure however that your reference sequence has a proper identifier in the fasta header as this will be used during analysis.

Nucleotide ambiguity codes in your reference fasta are supported but generally discouraged as this can have a negative effect on the primer removal process. You will get a warning during pre-processing if ambiguity nucleotides were found in your reference sequence.  

If you provide a reference fasta containing multiple sequences then analysis will be performed for every individual sequence in your input reference fasta, this will be split based on the reference identifier.  
This can be used in analysis of (for example) segmented viruses such as Influenza.


## Formatting your input GFF

Formatting of your input GFF is usually optional and very minimal.  

If you wish to have proper extraction and translation of aminoacids from the generated consensus sequence, then please make sure the input GFF has either a "Name" or an "ID" value in the attributes column of your input GFF.

For example: `Name="Nucleocapsid"` or `ID="N"` in the attributes column of your input GFF.

This information will be used to group the extracted aminoacid sequences of samples together during analysis.