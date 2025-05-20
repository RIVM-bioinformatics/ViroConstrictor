# Explanation of Amplicon Types

The term "amplicon type" describes a PCR amplicon and specifies which part of that amplicon is actually sequenced. This information is crucial because primers will be removed from the reads during analysis.

In the ViroConstrictor pipeline, primer removal is performed using our tool [AmpliGone](https://rivm-bioinformatics.github.io/AmpliGone/).

The information about the "amplicon type" is required because it tells the pipeline—and therefore AmpliGone—from which side(s) of a read the primer sequence should be removed.

Currently, there are three options: **"End-to-End"**, **"End-to-Mid."**, and **"Fragmented"**.

In short, **"End-to-End"** means that the sequenced read covers the full length of an amplicon, with the primer sequence present at both ends. In contrast, **"End-to-Mid"** means that the sequenced read only partially covers the length of the amplicon, so the primer sequence is present at only one end. **"Fragmented"** means that the amplicon is sequenced in multiple smaller reads, which may or may not contain primer sequences.

For more details, please refer to the [AmpliGone documentation regarding amplicon types](https://rivm-bioinformatics.github.io/AmpliGone/latest/amplicon-types/). Also, ensure that your sequencing and laboratory setup is configured for the correct amplicon type to achieve the best analysis results.
