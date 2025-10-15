# Multi-Reference Analysis

ViroConstrictor supports multiple reference sequences in a single analysis run, enabling flexible analysis workflows for complex datasets.

## When to Use Multi-Reference Analysis

Use multi-reference analysis when you have:
- Samples containing multiple pathogens
- Multiple strains of the same pathogen
- Uncertainty about which reference best fits your samples
- Segmented viruses with multiple subtypes (e.g., Influenza)

## Reference File Preparation

Multiple references must be combined into a single reference FASTA file containing all sequences you wish to use. This file can be provided via the `--reference` flag or through the samplesheet.

<!-- TODO: Add guidance on optimal number of references and performance considerations -->

---

## Analysis Modes

### 1. Full Analysis with Multiple References

**Use case**: Analyze samples against all provided references

**How it works**: ViroConstrictor performs a complete analysis for each reference in your FASTA file, generating separate outputs for each reference.

**Command**: Simply provide your multi-reference FASTA file using the `--reference` flag.

### 2. Best Reference Selection (Match-Ref)

**Use case**: Let ViroConstrictor automatically choose the best-fitting reference for each sample

**How it works**: ViroConstrictor evaluates all provided references and selects the one that best matches each sample's data.

**Requirements**:
- Multi-reference FASTA file
- Enable with `--match-ref` flag or set `Match-ref=True` in samplesheet

<!-- TODO: Explain the criteria used for "best" reference selection -->

!!! example "Multi-Reference FASTA Example"
    Example reference file with four measles virus subtypes for use with `--match-ref`:

    ```fasta
    >measles_subtype1  
    atgaaagtaaaactactggtcctg...  
    >measles_subtype2  
    atgagtcttctaaccgaggtcgaa...  
    >measles_subtype3  
    atgaacccaaatcaaaagataata...  
    >measles_subtype4  
    atgagtgacatcgaagccatggcg...
    ```

### 3. Segmented Virus Analysis

**Use case**: Analyze segmented viruses (like Influenza) where different segments may match different subtypes

**How it works**: ViroConstrictor selects the best-matching reference for each genome segment individually.

**Requirements**:
- Multi-reference FASTA with special header formatting
- Enable with both `--match-ref` and `--segmented` flags
- Or set both `Match-ref=True` and `Segmented=True` in samplesheet

#### Header Formatting for Segmented Viruses

**Required format**: `>{personal_identifier} {segment_name}|{segment_subtype}|{extra_information}`

**Key points**:
- `{segment_name}` must be consistent across all variants of the same segment
- This ensures proper result folder organization
- Headers are modified during analysis for correct output structure

<!-- TODO: Add validation tips for segmented virus header formatting -->


!!! example "Segmented Virus Reference Example"
    Example reference file for Influenza with three segments and three subtypes:

    ```fasta
    >A.HA_01 HA|H1|H1N1  
    atgaaagtaaaactactggtcc...  
    >A.HA_02 HA|H3|H3N2  
    atgaagactatcattgctttga...  
    >A.HA_03 HA|H5|H5N1  
    atgaagactatcattgctttga...  

    >A.MP_01 MP|MP|H1N1  
    atgagtcttctaaccgaggtcg...  
    >A.MP_02 MP|MP|H3N2  
    atgagccttcttaccgaggtcg...  
    >A.MP_03 MP|MP|H5N1  
    atgagtcttctaaccgaggtcg...  

    >A.NA_01 NA|N1|H1N1  
    atgaacccaaatcaaaagataa...  
    >A.NA_02 NA|N2|H3N2  
    atgaatccaaatcaaaagataa...  
    >A.NA_03 NA|N1|H5N1  
    atgaatccaaatcaaaagataa...
    ```

---

## Command Examples

### Basic Multi-Reference Analysis
```bash
viroconstrictor \
    --input /path/to/FASTQ-files \
    --output /path/to/output \
    --reference /path/to/multi-reference.fasta \
    --primers /path/to/primers.fasta \
    --features /path/to/features.gff \
    --platform nanopore \
    --amplicon-type end-to-end \
    --target viral-target
```

### Best Reference Selection
```bash
viroconstrictor \
    --input /path/to/FASTQ-files \
    --output /path/to/output \
    --reference /path/to/multi-reference.fasta \
    --primers /path/to/primers.fasta \
    --features /path/to/features.gff \
    --platform nanopore \
    --amplicon-type end-to-end \
    --target viral-target \
    --match-ref
```

### Segmented Virus Analysis
```bash
viroconstrictor \
    --input /path/to/FASTQ-files \
    --output /path/to/output \
    --reference /path/to/segmented-references.fasta \
    --primers /path/to/primers.fasta \
    --features /path/to/features.gff \
    --platform nanopore \
    --amplicon-type end-to-end \
    --target Influenza_A \
    --match-ref \
    --segmented
```

<!-- TODO: Add examples of samplesheet-based multi-reference analysis -->
