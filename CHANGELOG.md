# Changelog

## [1.0.0](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/compare/v0.4.2...v1.0.0) (2022-04-14)


### Bug Fixes

* Add option to not give primers on cli ([e4c054c](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/e4c054c570f670db6ce47bbfea788c4f937cc9ab))
* Check combination of samplesheet and cli args ([e4c054c](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/e4c054c570f670db6ce47bbfea788c4f937cc9ab))


### Code Refactoring

* move version to __init__ file ([e4c054c](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/e4c054c570f670db6ce47bbfea788c4f937cc9ab))


### Dependencies

* update AmpliGone to version 1.0.1 ([e4c054c](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/e4c054c570f670db6ce47bbfea788c4f937cc9ab))
* update TrueConsense version to v0.4.0 ([e4c054c](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/e4c054c570f670db6ce47bbfea788c4f937cc9ab))
* use Pypi for AmpliGone installation instead of building from source ([e4c054c](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/e4c054c570f670db6ce47bbfea788c4f937cc9ab))


### Documentation

* Add information for primer bed support ([e4c054c](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/e4c054c570f670db6ce47bbfea788c4f937cc9ab))
* update docs for multi-target analysis ([e4c054c](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/e4c054c570f670db6ce47bbfea788c4f937cc9ab))
* update documentation for new ViroConstrictor functions ([e4c054c](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/e4c054c570f670db6ce47bbfea788c4f937cc9ab))
* update documentation index page with correct required snakemake version ([e4c054c](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/e4c054c570f670db6ce47bbfea788c4f937cc9ab))
* update installation command for python3.10 ([e4c054c](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/e4c054c570f670db6ce47bbfea788c4f937cc9ab))

### [0.4.2](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/compare/v0.4.1...v0.4.2) (2022-04-12)


### Bug Fixes

* limit the self-updater to only update minor and patch versions ([f864f55](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/f864f55ab65e45c297a52e0f1cb9fa295c34a08f))

### [0.4.1](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/compare/v0.4.0...v0.4.1) (2022-03-08)


### Bug Fixes

* correct the memory assignment and threads for RemovePrimers job when no primers are given ([a9ebe02](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/a9ebe02cdc7995219d1a3c06fa9be6020881e43a))

## [0.4.0](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/compare/v0.3.0...v0.4.0) (2022-02-04)


### Features

* add an adjustable mismatch rate for the primer search ([7fc4265](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/7fc42651ac9bf9e926a93528a92cb9d6f58a9448))
* dynamically set usable cpus/threads in local execution mode ([1ec85bb](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/1ec85bb79c600d96c70bd4755451b174f25b0e80))
* Write a short report with the used/configured settings of the analysis ([0c46366](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/0c46366f72a8d1cd1a1ff87ad5e89234a59be9a5))


### Bug Fixes

* add missing --skip-updates flag ([f360161](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/f3601613b4bc007727239e03af5bf74feff8be65))
* corrected update function parameters ([7f37d90](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/7f37d90e310fde5767c3ef808f7ee5290d0fe6f6))
* finish an analysis when empty input fastq files are given ([d178587](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/d1785875d6e2b75f452003011d9b93554c69ccd1))
* properly exit with right exit codes and print closing message ([f1bbe13](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/f1bbe130c2ab00cad9d562a4e01113a4e33bc753))


### Performance Improvements

* tweak read qc filtering settings ([b3aab5a](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/b3aab5a095881668bf165be8c1eed0fb9a9c8a17))
* use adjusted mm2 mapping settings in nanopore-mode ([c027daa](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/c027daa31b37edff09ad285b9f0659b62a0328a7))


### Dependencies

* add FPDF version 1.7.2 to dependency lists ([cfcd9ee](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/cfcd9ee2863e923923cacbee19209f2f79fa6d43))
* add lenient version requirements to mamba ([07ebfcc](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/07ebfcc5cadf930a16f4f910ca990458e55043fd))
* add urllib3 v1.26 to dependency list for the auto-updater ([cfcd9ee](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/cfcd9ee2863e923923cacbee19209f2f79fa6d43))
* change biopython version in python setupfile ([07ebfcc](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/07ebfcc5cadf930a16f4f910ca990458e55043fd))
* change pysam to lenient versioning ([cfcd9ee](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/cfcd9ee2863e923923cacbee19209f2f79fa6d43))
* inlcude conda (lenient) version 4.11.x in main env for compatibility ([cfcd9ee](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/cfcd9ee2863e923923cacbee19209f2f79fa6d43))
* pin pysamstats version to 1.1.2 ([cfcd9ee](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/cfcd9ee2863e923923cacbee19209f2f79fa6d43))
* pin tqdm (lenient) version to 4.62.x ([cfcd9ee](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/cfcd9ee2863e923923cacbee19209f2f79fa6d43))
* pin TrueConsense to version 0.3.0 ([cfcd9ee](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/cfcd9ee2863e923923cacbee19209f2f79fa6d43))
* remove pysamstats from 'clean' environment ([cfcd9ee](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/cfcd9ee2863e923923cacbee19209f2f79fa6d43))
* update AmpliGone 0.4.0 --> 0.4.1 ([ea06728](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/ea06728c3e2c05e5adba4320683d75c59a1d3a16))
* update AmpliGone 0.4.1 --> 0.4.3 ([af62d22](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/af62d2294d48249a8108ce03bdaae72f626cdb51))
* update AmpliGone version 0.2.3 --> 0.4.0 ([cfcd9ee](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/cfcd9ee2863e923923cacbee19209f2f79fa6d43))
* update bcftools 1.12 --> 1.14 ([cfcd9ee](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/cfcd9ee2863e923923cacbee19209f2f79fa6d43))
* update bedtools 2.29.2 --> 2.30.0 ([cfcd9ee](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/cfcd9ee2863e923923cacbee19209f2f79fa6d43))
* update biopython 1.78 --> 1.79 ([cfcd9ee](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/cfcd9ee2863e923923cacbee19209f2f79fa6d43))
* update fastp 0.20.1 --> 0.23.2 ([cfcd9ee](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/cfcd9ee2863e923923cacbee19209f2f79fa6d43))
* update fastqc 0.11.8 --> 0.11.9 ([cfcd9ee](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/cfcd9ee2863e923923cacbee19209f2f79fa6d43))
* update minimap2 & mappy 2.17 --> 2.24 ([cfcd9ee](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/cfcd9ee2863e923923cacbee19209f2f79fa6d43))
* update multiqc 1.9 --> 1.11 ([cfcd9ee](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/cfcd9ee2863e923923cacbee19209f2f79fa6d43))
* update pandas 1.2.3 --> 1.3.5 ([cfcd9ee](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/cfcd9ee2863e923923cacbee19209f2f79fa6d43))
* update parmap 1.5.2 --> 1.5.3 ([cfcd9ee](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/cfcd9ee2863e923923cacbee19209f2f79fa6d43))
* update PyYaml 5.4.1 --> 6.0 ([cfcd9ee](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/cfcd9ee2863e923923cacbee19209f2f79fa6d43))
* update samtools 1.11 --> 1.14 ([cfcd9ee](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/cfcd9ee2863e923923cacbee19209f2f79fa6d43))
* update seqkit 0.14.0 --> 2.1.0 ([cfcd9ee](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/cfcd9ee2863e923923cacbee19209f2f79fa6d43))
* update snakemake 6.4.1 --> 6.13.1 ([cfcd9ee](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/cfcd9ee2863e923923cacbee19209f2f79fa6d43))


### Documentation

* add importance of primer name and number to docs ([7107398](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/7107398833b1082b4d5c188bf6edbed47ee73996))
* include primer mismatch rate flag to documentation ([531f526](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/531f5262dcad16ad9386c59992dc2f64d3f869c9))
* update installation instructions ([5944d22](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/5944d225361797ef2c12b8cb2e5dc3c0cbd12c85))

## [0.3.0](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/compare/v0.2.3...v0.3.0) (2021-11-15)


### Features

* allow ViroConstrictor to update itself to latest release ([11ad616](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/11ad6160f7748556598e7447b35975b9dfb378fd))
* flexible memory requirements per workflow job ([7c189cd](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/7c189cdac758f49fcc578ec228abb21c2f5af1ad))


### Bug Fixes

* change consensus calling to a 'mid' cpu-task ([3094b78](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/3094b78bbce88f9df5b71e8de8ac2c5ce737374a))

### [0.2.3](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/compare/v0.2.2...v0.2.3) (2021-11-03)


### Documentation

* add explanatory comments to docs config ([f3e05ec](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/f3e05ecbc47e7f68edefa8829ce95876b147241e))

### [0.2.2](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/compare/v0.2.1...v0.2.2) (2021-11-03)


### Bug Fixes

* change amplicon covs script to complete with empty input [#5](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/issues/5) ([6900638](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/6900638375ad95b96771de15bbb2ea7fa9d3ba14))
* complete workflow with primers set to NONE ([47becf8](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/47becf85eba7c92f70df27fe93ccffdf27efae1c))
* use 'map-ont' mm2 preset in nanopore workflow ([f89d5cf](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/f89d5cf7154c139eab1ce361fc4107c4cd568217))
