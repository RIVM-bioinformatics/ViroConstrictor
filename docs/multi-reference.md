# Working with multiple references

ViroConstrictor is able to work with multiple references in a single analysis run. This can be beneficial when you're working with samples that contain multiple pathogens, multiple strains of the same pathogen, or when you're not sure which reference to use for a specific sample.

Multiple references are provided in a singular reference fasta file. This file should contain all the references you want to use in the analysis. The reference fasta file should be provided to ViroConstrictor with the `--reference` flag or through the samplesheet.


## Run a full analysis with multiple references

When you want to run a full analysis with multiple references you can provide the reference fasta file with the `--reference` flag. The reference fasta file should contain all the references you want to use in the analysis.

By default, the full analysis will now be ran for each reference in the provided reference fasta file. This will result in multiple results for each provided reference.

## Choose the best reference for each sample

If you have a sequencing protocol is effective for multiple strains of the same pathogen, and you're not sure which reference to use for a specific sample, then you can let ViroConstrictor choose the best reference for each sample out of a larger set of references.

To do this provide all potential references in a single reference fasta file through the `--reference` flag or through the samplesheet. Additionally, provide the `--match-ref` flag or set the `Match-ref` column in the samplesheet to `True`.

When the `--match-ref` flag is provided or the `Match-ref` column is set to `True` in the samplesheet, ViroConstrictor will try to choose the best reference for each sample out of the provided references.

!!! example "Example of a reference fasta file"
    Below is an example of a reference fasta file depicting 4 subtypes of the measles virus. This file can be provided to ViroConstrictor with the `--reference` flag or through the samplesheet. With `--match-ref` enabled or the `Match-ref` column set to `True` in the samplesheet, ViroConstrictor will pick the best matching reference out of these 4 subtypes, this will be done individually for each sample.

    >   \>measles_subtype1  
        atgaaagtaaaactactggtcctg...  
        \>measles_subtype2  
        atgagtcttctaaccgaggtcgaa...  
        \>measles_subtype3  
        atgaacccaaatcaaaagataata...  
        \>measles_subtype4  
        atgagtgacatcgaagccatggcg...

## Choose the best reference for each sample with segmented viruses

If you have a sequencing protocol that is effective for multiple strains of the same virus, and that virus has a segmented genomic structure, then you can let ViroConstrictor choose the best fitting reference for each segment of the virus for each sample. This can be beneficial when you're working with segmented viruses like influenza.

This can be achieved by providing a single reference fasta file with all the reference segments, and variations of the segments if applicable, of the virus. To let ViroConstrictor choose the best fitting reference for each segment of the virus for each sample, provide both the `--match-ref` flag and the `--segmented` flag. Or set both the `Match-ref` and `Segmented` columns in the samplesheet to `True`.

Additionally, extra formatting of the reference fasta file is required. Each fasta-header should follow a format as shown below:

>\>{personal_identifier} {segment_name}|{segment_subtype}|{extra information}

!!! example "Example of a reference fasta file for segmented viruses"
    Below is an example of a reference fasta file depicting 3 segments of 3 subtypes of the influenza virus.  
    This file can be provided to ViroConstrictor with the `--reference` flag or through the samplesheet. With `--match-ref` and `--segmented` enabled or the `Match-ref` and `Segmented` columns set to `True` in the samplesheet, ViroConstrictor will pick the best matching reference for each segment of the virus for each sample.

    >   \>A.HA_01 HA|H1|H1N1  
    >   atgaaagtaaaactactggtcc...  
    >   \>A.HA_02 HA|H3|H3N2  
    >   atgaagactatcattgctttga...  
    >   \>A.HA_03 HA|H5|H5N1  
    >   atgaagactatcattgctttga...  

    >   \>A.MP_01 MP|MP|H1N1  
    >   atgagtcttctaaccgaggtcg...  
    >   \>A.MP_02 MP|MP|H3N2  
    >   atgagccttcttaccgaggtcg...  
    >   \>A.MP_03 MP|MP|H5N1  
    >   atgagtcttctaaccgaggtcg...  
    >    
    >   \>A.NA_01 NA|N1|H1N1  
    >   atgaacccaaatcaaaagataa...  
    >   \>A.NA_02 NA|N2|H3N2  
    >   atgaatccaaatcaaaagataa...  
    >   \>A.NA_03 NA|N1|H5N1  
    >   atgaatccaaatcaaaagataa...  

    Please note that after choosing the best fitting reference for each segment of the virus for each sample, the fasta-header will be modified to ensure that the results folder is structured correctly. This requires that the {segment-name} will be the same for every segment-variant provided in the reference fasta file.