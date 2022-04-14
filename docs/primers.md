# Formatting your input primers

In order to get optimal results, please make sure the fasta headers in your fasta file with primers are formatted properly. Please make sure the fasta headers are formatted according the following format:

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