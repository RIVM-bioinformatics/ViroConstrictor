# Amplicon Types

The "amplicon type" describes the relationship between PCR amplicons and how they are sequenced. This information is crucial for proper primer removal during analysis, as it determines from which side(s) of a read the primer sequences should be removed.

ViroConstrictor uses [AmpliGone](https://rivm-bioinformatics.github.io/AmpliGone/) for primer removal, which requires accurate amplicon type specification to function correctly.

## Available Amplicon Types

ViroConstrictor supports three amplicon types: **End-to-End**, **End-to-Mid**, and **Fragmented**. The choice depends on your sequencing protocol and how reads cover the amplicons.

---

## End-to-End

**Definition**: Reads cover the full length of an amplicon from primer to primer.

**Characteristics**:

* Primers are present at both ends of reads
* Reads span the complete amplicon length
* Both forward and reverse primers need removal

**Common platforms**:

* Nanopore sequencing
* Illumina MiSeq (long amplicons)

**Primer removal logic**:

* Both forward and reverse reads should start at a forward primer and end at a reverse primer
* Read orientation doesn't affect primer detection
* Primers removed from both ends of each read

<!-- TODO: Add guidance on minimum/maximum amplicon lengths for end-to-end -->
<!-- TODO: Add explanatory images to visualize how this works -->

---

## End-to-Mid

**Definition**: Reads partially cover amplicons and overlap at the amplicon midpoint.

**Characteristics**:

* Primers present at only one end of reads
* Forward reads start at forward primers
* Reverse reads start at reverse primers
* Reads meet/overlap in the middle of amplicons

**Common platforms**:

* Illumina MiSeq (standard paired-end)
* Most paired-end short-read platforms

**Primer removal logic**:

* Forward reads: primer removal from 5' end only
* Reverse reads: primer removal from 5' end only (reverse primer)
* No primer removal from 3' ends

<!-- TODO: Explain how overlap quality affects consensus calling -->
<!-- TODO: Add explanatory images to visualize how this works -->

---

## Fragmented

**Definition**: Multiple smaller reads cover single larger amplicons, with reads potentially scattered across the amplicon length.

**Characteristics**:

* Reads may or may not contain primer sequences
* Not all reads start at primer positions
* Reads can be internal to amplicons
* Requires proximity-based primer detection

**Common platforms**:

* Illumina NextSeq (shorter reads, longer amplicons)
* Short-read platforms with large amplicons

**Primer removal logic**:

* Primers removed only from reads directly linked to primer regions
* Uses "fragment lookaround size" (default: 10bp) to determine primer proximity
* Reads starting/ending within 10bp of primer positions are considered primer-containing

**Advanced option**: The fragment lookaround size can be adjusted with `--fragment-lookaround-size` flag if needed.

<!-- TODO: Add guidance on optimal fragment lookaround size for different protocols -->

---

## Choosing the Right Amplicon Type

### Decision Guide

1. **Do your reads span entire amplicons?** → End-to-End
2. **Do paired reads meet in the middle of amplicons?** → End-to-Mid  
3. **Do you have short reads covering long amplicons?** → Fragmented

### Platform Guidelines

| Platform | Typical Amplicon Type | Notes |
|----------|----------------------|-------|
| Nanopore | End-to-End | Long reads usually span full amplicons |
| Illumina MiSeq | End-to-Mid | Standard paired-end overlap protocol |
| Illumina NextSeq | Fragmented | Shorter reads may require fragmented approach |

!!! warning "Protocol Matters More Than Platform"
    While these are typical associations, your specific sequencing protocol and amplicon design determine the correct type. A MiSeq could use End-to-End for very short amplicons, or NextSeq could use End-to-Mid with appropriate amplicon sizes.

---

## Why Not Use Standard Terminology?

Traditional terms like "paired-end/single-end" or "short/long read" don't always convey the specific information AmpliGone needs for primer removal. 

**Example scenarios**:

* Short reads might have primers on both ends (End-to-End) or just one end (End-to-Mid)
* Platform alone doesn't determine amplicon coverage patterns
* Experimental design varies independently of sequencing technology

The amplicon type system provides precise information about primer locations without requiring multiple command-line arguments or platform-specific assumptions.

<!-- TODO: Add troubleshooting section for common amplicon type misassignment issues -->
