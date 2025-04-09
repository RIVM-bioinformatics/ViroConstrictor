# Working with Multiple References

ViroConstrictor can work with multiple references in a single analysis run. This is useful when you are working with samples that contain multiple pathogens, multiple strains of the same pathogen, or when you are unsure which reference to use for a specific sample.

Multiple references should be combined into a single reference FASTA file. This file must contain all the references you wish to use in the analysis. The reference FASTA file can be provided to ViroConstrictor either with the `--reference` flag or via the samplesheet.

## Run a Full Analysis with Multiple References

To perform a full analysis with multiple references, simply provide the reference FASTA file using the `--reference` flag. The file should include all the references you want to use.  
By default, the full analysis will be run for each reference in the provided FASTA file, resulting in multiple outputs—one for each reference.

## Choose the Best Reference for Each Sample

If your sequencing protocol is effective for multiple strains of the same pathogen and you are unsure which reference to use for a particular sample, you can have ViroConstrictor choose the best reference for each sample from a larger set of references.

To enable this, provide all potential references in a single reference FASTA file via the `--reference` flag or in the samplesheet. Additionally, supply the `--match-ref` flag or set the `Match-ref` column in the samplesheet to `True`.  
When enabled, ViroConstrictor will attempt to select the best reference for each sample from those provided.

!!! example "Example of a Reference FASTA File"
    Below is an example of a reference FASTA file depicting four subtypes of the measles virus. This file can be provided to ViroConstrictor via the `--reference` flag or through the samplesheet. With `--match-ref` enabled or the `Match-ref` column set to `True` in the samplesheet, ViroConstrictor will choose the best matching reference for each sample.

    >   \>measles_subtype1  
        atgaaagtaaaactactggtcctg...  
    >   \>measles_subtype2  
        atgagtcttctaaccgaggtcgaa...  
    >   \>measles_subtype3  
        atgaacccaaatcaaaagataata...  
    >   \>measles_subtype4  
        atgagtgacatcgaagccatggcg...

## Choose the Best Reference for Each Sample with Segmented Viruses

If you have a sequencing protocol that works for multiple strains of the same virus and that virus has a segmented genome (such as influenza), you can have ViroConstrictor choose the best fitting reference for each segment for each sample.

This is achieved by providing a single reference FASTA file that contains all the reference segments (and their variants, if applicable) of the virus. To enable the selection of the best reference for each segment, provide both the `--match-ref` flag and the `--segmented` flag, or set both the `Match-ref` and `Segmented` columns in the samplesheet to `True`.

Additionally, extra formatting of the reference FASTA file is required. Each FASTA header should follow the format shown below:

>\>{personal_identifier} {segment_name}|{segment_subtype}|{extra information}


!!! example "Example of a Reference FASTA File for Segmented Viruses"
    Below is an example of a reference FASTA file depicting three segments of three subtypes of the influenza virus.  
    This file can be provided to ViroConstrictor via the `--reference` flag or through the samplesheet. With both `--match-ref` and `--segmented` enabled (or the corresponding columns set to `True` in the samplesheet), ViroConstrictor will select the best matching reference for each segment for each sample.

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

    >   \>A.NA_01 NA|N1|H1N1  
    >   atgaacccaaatcaaaagataa...  
    >   \>A.NA_02 NA|N2|H3N2  
    >   atgaatccaaatcaaaagataa...  
    >   \>A.NA_03 NA|N1|H5N1  
    >   atgaatccaaatcaaaagataa...

Please note that after choosing the best fitting reference for each segment of the virus for each sample, the FASTA header will be modified to ensure that the results folder is structured correctly. This requires that the `{segment_name}` format is consistent for every segment variant provided in the reference FASTA file.
