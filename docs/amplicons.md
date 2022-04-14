# Explanation of Amplicon Types

The "amplicon type" is basically a description of a PCR amplicon, and which part of said amplicon is actually sequenced.  
This information is important because given primers will be removed from the reads during the analysis.

In the ViroConstrictor pipeline the primer-removal is done by our tool [AmpliGone](https://rivm-bioinformatics.github.io/AmpliGone/).

The information regarding the "amplicon type" is required because it tells te pipeline, and therefore also AmpliGone, from which side(s) of a read the primer-sequence has to be removed.

Currently there are two options: "End-to-End" or "End-to-Mid".  

In short: "End-to-End" means that the sequenced read covers the full length of an amplicon. Meaning that the primer-sequence is present at both ends of a read.  
"End-to-Mid" means that the sequences read only *partially* covers the length of the amplicon. Meaning that the primer-sequence is present at only one end of a read.

Also see the [AmpliGone documentation regarding amplicon-types](https://rivm-bioinformatics.github.io/AmpliGone/latest/amplicon-types/) where we go into more detail regarding these types and what these terms mean.  
Please check your sequencing and laboratory setup for the amplicon type to match your analysis to ensure the best results.