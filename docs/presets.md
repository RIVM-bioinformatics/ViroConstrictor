# Analysis Presets

ViroConstrictor includes a preset system (available since [version 1.2.0](changelog.md#120-2023-01-18)) that provides optimized analysis settings for specific pathogens.

## What Are Presets?

Presets are predefined analysis configurations tailored to specific viral families or pathogens. Different viruses require different analysis approaches - for example, segmented viruses like Influenza need different settings than non-segmented viruses like SARS-CoV-2.

## How Presets Work

### Automatic Selection
- **Input matching**: ViroConstrictor matches your viral target name to the closest preset
- **Fuzzy matching**: Exact matches aren't required - the system finds the best fit
- **Confidence threshold**: Matches below 40% confidence use DEFAULT settings
- **Fallback**: Unknown targets automatically use DEFAULT preset

### Usage Options
- **Command-line**: Use `--target` or `--preset` flag for entire analysis
- **Samplesheet**: Specify individual presets via the "Virus" column
- **Disable**: Use `--disable-presets` / `-dp` to force DEFAULT settings

<!-- TODO: Add examples of successful and failed preset matching -->

!!! info "When to Disable Presets"
    Use `--disable-presets` / `-dp` to force DEFAULT settings when:
    
    - Analyzing experimental or in-development sequencing protocols
    - Working with pathogens not covered by existing presets
    - Troubleshooting analysis issues
    - Requiring consistent settings across diverse viral targets

---

## Available Presets

!!! warning "Version-Dependent Presets"
    Preset availability and settings may change between ViroConstrictor versions. The information below corresponds to the current documentation version shown at the top of this page.

### Preset Matching Examples

ViroConstrictor uses fuzzy matching, so exact names aren't required:

* `--target influenza_a_h3n2` → matches **INFLUENZA** preset via `INFLUENZA_A` alias.
* `--target sars-cov-2` → matches **SARSCOV2** preset via `SARS-COV-2` alias.
* `--target measles_outbreak` → matches **PARAMYXOVIRIDAE** preset via `MEASLES` alias.

<!-- TODO: Add more examples of edge cases and unexpected matches --> 

### Preset Details

| Preset                | Preset Aliases for matching  | Optimizations |
|-----------------------|------------------------------|---------------|
| **DEFAULT**           | N/A                          | Standard ViroConstrictor settings for general viral analysis.<br>Used when no specific preset matches or when presets are disabled. |
| **SARSCOV2**          | SARS-COV-2, SARS2, COVID, COV, SARSCOV, CORONAVIRUS | Optimized read-alignment settings for SARS-CoV-2 reference mapping.<br>Tailored for SARS coronavirus 2 genome analysis. |
| **INFLUENZA**         | FLU, INF, INFLU, INFLUENZA, INFLUENZA_A, INFLUENZA_B | Enhanced read-filtering and alignment-filtering for segmented Influenza genomes.<br>Built-in compatibility for MBTuni-12/13 universal influenza primer amplicons.[^1] |
| **PARAMYXOVIRIDAE**   | MEASLES, MUMPS, MEV, MUV, MEASLES_VIRUS, MUMPS_VIRUS, PARAMYXOVIRUS, PARAMYXOVIRIDAE, MORBILLIVIRUS, RUBULAVIRINAE, ORTHORUBULAVIRUS | Specialized alignment-filtering for Paramyxoviridae family viruses with optimized settings to cover the non-coding region between Matrix and Fusion genes.<br>Covers Measles, Mumps, and related viruses. |
| **HEPATOVIRUS**       | HEPATOVIRUS, HEPATITIS_A, HEPATITIS_A_VIRUS, HAV | Alignment-filtering optimized for Hepatovirus A. Includes specialized alignment-filtering settings for improved genome recovery to compensate for sequencing artefacts in Picornaviridae sequencing data.<br>Currently only tailored for Hepatovirus A. |
| **PNEUMOVIRIDAE**     | PNEUMOVIRIDAE, RSV, RESPIRATORY_SYNCYTIAL_VIRUS, RESPIRATORY_SYNCYTIAL_VIRUS_A, RESPIRATORY_SYNCYTIAL_VIRUS_B, HMPV, METAPNEUMOVIRUS, HUMAN_METAPNEUMOVIRUS, PNEUMOVIRUS | Optimized for respiratory syncytial virus (RSV) and human metapneumovirus (HMPV) analysis.<br>Enhanced settings for spliced genome alignment and virtual primer support. |
| **ORTHOHEPADNAVIRUS** | HBV, HEPB, HEPATITIS_B, HEPATITIS_B_VIRUS, ORTHOHEPADNAVIRUS | Specialized settings for Hepatitis B virus analysis.<br>Configured with increased minimum read length requirements for improved genome recovery. |
| **ENTEROVIRUS**       | ENTEROVIRUS, EV, ECHO, ECHOVIRUS, COX, COXSACKIE, COXSACKIEVIRUS, RHINOVIRUS, RHINO, POLIO, POLIOVIRUS | Enhanced read-filtering and alignment settings for Enterovirus family analysis.<br>Covers Poliovirus, Coxsackievirus, Echovirus, and Rhinovirus. |

---

## Important Considerations

### Matching Warnings

ViroConstrictor displays warnings when:
- **Low confidence matches**: Your input has <40% similarity to any preset
- **Unexpected matches**: Your target matches a potentially unintended preset
- **No match found**: Your input doesn't match any preset sufficiently

!!! tip "Best Practices"
    - Run with `--dryrun` first to verify preset assignments
    - Check warnings carefully after analysis completion
    - Use `--disable-presets` if matching seems incorrect

<!-- TODO: Add guidance on interpreting preset matching confidence scores -->


### Requesting New Presets

Need a preset for a pathogen not currently supported? [Submit a feature request on our GitHub](https://github.com/RIVM-bioinformatics/ViroConstrictor/issues/new/choose) with details about:

- Target pathogen or viral family
- Sequencing protocol used
- Specific analysis challenges
- Expected benefits from a custom preset

---

[^1]: [Zhou B, Donnelly ME, Scholes DT, et al. Single-Reaction Genomic Amplification Accelerates Sequencing and Vaccine Production for Classical and Swine Origin Human Influenza A Viruses. Journal of Virology. 83, 10309-10313.](https://doi.org/10.1128/JVI.01109-09)