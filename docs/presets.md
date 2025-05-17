# Working with Presets

Since [version 1.2.0](changelog.md#120-2023-01-18) ViroConstrictor has a 'preset' system in place.
The primary goal of this system is to provide tailored analysis settings and methods specific to the pathogen or group of pathogens being analyzed.
As an example, even though Influenza and SARS-CoV-2 are both respiratory viruses, they require some different analysis settings to ensure the best possible results as the viral genomes and the sequencing protocols are very different (segmented vs non-segmented).

A preset can be set both for an entire analysis run with the `--target/--preset` flag or on a per sample basis through the "Virus" column in the samplesheet.
The input value does not have to completely correspond to a defined preset, the closest matching preset is chosen unless there is not enough certainty (<40%) to match the input to a preset.
If there is not enough certainty to match the input to a preset then the "DEFAULT" preset/settings will be used as a fallback.

!!! info "Disabling the use of presets"
    The use of presets can be disabled by providing the `--disable-presets/-dp` flag. This will make sure the "DEFAULT" preset/settings will be used for all samples given in a run.  
    Please refer to the [full list of commandline options](manual.md#overview-of-command-line-options) for more information.

    Disabling the presets is especially useful if you're trying to analyze an in-progress sequencing protocol, or if you're trying to analyze data that cannot be associated with an existing preset [listed in the table below](presets.md#currently-available-presets)

## Currently available presets

Below you can find a brief summary of the currently available presets within ViroConstrictor.  
!!! warning "Presets can change throughout ViroConstrictor versions"
    Please note that additions, removals or changes in presets happen in ViroConstrictor releases.  
    These docs correspond to a version of ViroConstrictor, the version of ViroConstrictor/docs that you are viewing right now can be seen at the top of this page. 

| Preset    | Aliases | Notes |
|-----------|-------|---------|
| DEFAULT   | N/A | This is the default analysis mode for ViroConstrictor.<br>Read-filtering, read-alignment and consensus-calling will all use default settings |
| SARSCOV2  | <ul><li>SARS-COV-2</li><li>SARS2</li><li>COVID</li><li>COV</li><li>SARSCOV</li><li>CORONAVIRUS</li></ul> | This preset is specific for the analysis of SARS-CoV-2 data.<br>Adjustments have been made to read-alignment settings for optimal mapping of reads against the reference genome. |
| INFLUENZA | <ul><li>FLU</li><li>INF</li><li>INFLU</li><li>INFLUENZA</li><li>INFLUENZA_A</li><li>INFLUENZA_B</li></ul> | This preset is specific for the analysis of Influenza A and Influenza B data.<br>There are specific settings for read-filtering, and alignment-filtering to assure the best result for the segmented influenza genomes.<br>The settings are chosen with the amplicons produced by the MBTuni-12/13 universal influenza primers in mind.[^1] |

Presets are set based on the closest matching input-value relative to a preset alias, a 100% match is therefore not necessary.  
As an example when the input `--target influenza_a_h3n2` is given, it will be matched to the alias `INFLUENZA_A` which corresponds to the "INFLUENZA" preset.

If there's a significant distance between your given input and the matched preset then ViroConstrictor will display a warning for this. This warning is shown **after** the analysis (or dryrun) is completed to ensure the warning stays clearly readable and is not pushed off-screen.

!!! warning "Please beware of unintentional preset matchings"
    It may be possible that a preset gets assigned to a provided target unintentionally, as the preset matching happens based on your user-input, and because the closest matching alias will be used.  
    Please inspect the possible preset-related warning carefully to see if this is the case.

    We recommend running ViroConstrictor with the `--dryrun` flag first to make sure presets get assigned correctly.

    If any presets are not assigned correctly, consider running ViroConstrictor with the `--disable-presets` flag to always use the "DEFAULT" preset/settings.


If you wish for another pathogen to get its own preset, or if you're working with the analysis of a specific virus that would benefit from a preset, please [request this through an issue on our GitHub](https://github.com/RIVM-bioinformatics/ViroConstrictor/issues/new/choose).


[^1]: [Zhou B, Donnelly ME, Scholes DT, et al. Single-Reaction Genomic Amplification Accelerates Sequencing and Vaccine Production for Classical and Swine Origin Human Influenza A Viruses. Journal of Virology. 83, 10309-10313.](https://doi.org/10.1128/JVI.01109-09)