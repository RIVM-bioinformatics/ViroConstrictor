# Changelog

## [1.6.6](https://github.com/RIVM-bioinformatics/ViroConstrictor/compare/v1.6.5...v1.6.6) (2026-03-01)


### Bug Fixes

* change log.warn to log.warning for deprecationwarning ([c39f41e](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/c39f41e281daf259ebd0d2d5314121c4ad77cbf4))
* enhance primer name parsing to support insert size format ([4e4eff3](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/4e4eff3835f1e2851c533cc92d1ca7c1bc23b70f))
* handle missing required assets in CLIparser ([c39f41e](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/c39f41e281daf259ebd0d2d5314121c4ad77cbf4))

## [1.6.5](https://github.com/RIVM-bioinformatics/ViroConstrictor/compare/v1.6.4...v1.6.5) (2026-02-23)


### Bug Fixes

* add url parsing and audit before opening the url in updater process ([7f221d0](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/7f221d08241c6e4c9f356a81b5d42930ed583316))
* expand primer regex patterns to include artic v5.4.2 scheme ([2b898e8](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/2b898e8d9752af6d99ed8652938fd13f826194cb))
* improve the updater-process to use no longer use the mamba api in order to circumvent solver errors ([cb76f34](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/cb76f34730efd4c6f7737613417003af24f9ccd0))
* update primer regex pattern to allow for optional primer version in naming ([6b165e4](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/6b165e446f468813b0c5bdb2f6eb714167576582))
* update repo_root path to reference the correct parent directory ([7a98e8c](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/7a98e8ce793b358aa1c5300da8370b4854e3ef96))


### Documentation

* add docstrings in numpy style to `test_group_refs.py` ([a37d156](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/a37d156457844426e6b275fc9467d98639c42064))
* add docstrings to `test_update.py` in numpy docstring style ([7b86595](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/7b86595452ddcbe54ab3aaaca3e638e03f0bed4e))
* added and/or updated all docstrings in `update.py` in numpy docstrings style ([b59f7c8](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/b59f7c897ee2a98e337263e6dd10f9b05ca15b00))
* fix typo in presets documentation ([18c2ecc](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/18c2ecc60989744b8f3ca0dfea1574eb5cfd298b))
* update docstrings and formatting in update.py ([8b14e12](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/8b14e12d6abcdd68245bcbfca44686d2420de8a3))
* update links in manual.md to point to existing files ([33167d9](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/33167d97ca8f01f1bab194503a35ddb07f0a4840))

## [1.6.4](https://github.com/RIVM-bioinformatics/ViroConstrictor/compare/v1.6.3...v1.6.4) (2026-02-10)


### Features

* add a combined results folder with aggregates ([f52a2dc](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/f52a2dcf046462afb9e5307142dd75f57c3540e7))
* add combined results structure for more human friendly interpretation when dealing with many references. ([246e3ed](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/246e3ed79b4c46a0dccf787ae9b2fd0ff99d0c31))
* add enterovirus preset aliases and parameters for workflow(s) ([86c09a3](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/86c09a36dd07a834ca2418f19c3019f934755fd2))
* add functionality to combine amino acid results by sample and virus, including helper functions for feature extraction ([3e8507f](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/3e8507fe9d697ff8c6e5b0e73899b59c7a2c87c7))
* add hbv preset ([89089f7](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/89089f74d26a6bf4c31c490be1fb75d28d868839))


### Bug Fixes

* add various small changes in unit test and validation behaviour ([2b16a8b](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/2b16a8b0696822be9a005c1eaa738f4a0b2e5edc))
* adjust time allocations to accommodate for O(n^2) situation for (very) large samples ([d4418dd](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/d4418ddee18e596e5d62150368b44fe42ede2a08))
* amplicon_covs - zero-pad amplicon numbers to 3 digits for proper sorting ([197e339](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/197e3399341a9353ef814d3e5796502e21271049))
* enhance amplicon start and end calculations to ensure overlapping regions are excluded ([080722a](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/080722ae9dc23342952e9fa0db5d8c1138d088d1))
* ensure feature lists are constructed on a per-sample basis in `construct_rule_all` ([515ae64](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/515ae643403552504d4ff44367b3d77935cea72f))
* ensure mean coverage per amplicon is calculated excluding the amplicon-overlaps ([b2e4d9d](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/b2e4d9d129469a1a39e75272abf25affecbbc2c1))
* ensure NaN values in Primer_file and Feat_file columns are replaced with string "NONE" ([7352a64](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/7352a64c864b4d70b13785f788247313cc28de78))
* escape parameters in virus and refid lists for proper handling ([f4fefd9](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/f4fefd92c55c86ac2afb7723524b382700861517))
* refactor make_pickle rule to simplify pickle creation and add logging in concat_aminoacids rule ([5e35714](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/5e357144fe6bce408c04fec391a450e3c55aae16))
* remove quoting for virus and refid parameters in aggregation rules which had an unintentional side effect ([25f7674](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/25f76748d29ffb0368257d65da26ba0442f1c26d))
* resolve header row handling for mutations and width_of_coverage files ([036a2be](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/036a2be2c875b1b4aaff0dfc1bbc96e90e54d81f))
* round mean coverage calculation to two decimal places. ([fbf43f9](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/fbf43f993b060c72ed7e0ef4ba16416ac950ccdd))
* update `combine_aminoacids_by_sample` rule to handle features correctly ([5135558](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/513555861af67363ebb5d0090b8925946699b946))
* update amino acid feature checks to use correct variable references in aggregation rules ([c531e02](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/c531e0292a0fa319c815ea6145cc3c5c6812230b))
* update amplicon start and end calculations to ensure overlapping regions are excluded ([b38a8d0](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/b38a8d04c0470077a644702b0bc78dc3eb37e28d))
* update CSV reading to handle NA values correctly in AggregateCombinedFiles and CombineTabular ([a534def](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/a534def53092b0d4d6dde9b3de1ca2e4e126dd9c))
* update fasta header formatting to always include Virus and RefID information ([3073dea](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/3073deada1618b6a00b79727b9f90b461db865ae))
* update runtime function tests to reflect quadratic scaling ([c5d1da7](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/c5d1da742870e629322d4f008765949944f462a5))


### Documentation

* addition of docstrings in `aggregate_combined_files.py` ([6002df5](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/6002df565108d876bb448c839896f73172510131))
* addition of docstrings in `combine_fasta.py` ([f91be72](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/f91be722392148ddc13a55e060f157195ebcd98f))
* addition of docstrings in `combine_tabular.py` ([63b761b](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/63b761b50417635cd1038da9257ffb39ac2eac77))
* addition of docstrings in `extract_sample_from_fasta.py` ([62461a8](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/62461a8c30abefa76ade240ed1d150deb4f03ffd))
* update presets documentation with new virus categories and optimized settings ([2b7c4cd](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/2b7c4cdf8918ac6c4ef72bce922df183bf3adb43))


### Miscellaneous Chores

* **release:** empty commit ([9a0aeba](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/9a0aeba04d60e5eeb7c0123a85c1dafe20a3f45b))

## [1.6.3](https://github.com/RIVM-bioinformatics/ViroConstrictor/compare/v1.6.2...v1.6.3) (2025-12-29)


### Bug Fixes

* add additional primer name option (MuV-NGS_19-alt22_LEFT) ([#169](https://github.com/RIVM-bioinformatics/ViroConstrictor/issues/169)) ([94056f8](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/94056f87c338c028e27beafc0cfd9cb1aace3c23))

## [1.6.2](https://github.com/RIVM-bioinformatics/ViroConstrictor/compare/v1.6.1...v1.6.2) (2025-12-17)


### Bug Fixes

* small fixes related to samples df ([#166](https://github.com/RIVM-bioinformatics/ViroConstrictor/issues/166)) ([b018472](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/b018472ab443add8ff0584e9f820addeddf92a91))

## [1.6.1](https://github.com/RIVM-bioinformatics/ViroConstrictor/compare/v1.6.0...v1.6.1) (2025-12-12)


### Bug Fixes

* adjust the SARS-CoV-2 preset for proper recovery of +800bp deletion in ORFs 7a, 7b, and 8 ([437fc09](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/437fc09bb829d8fa67f333ef9550fbfaebf70bf9))

## [1.6.0](https://github.com/RIVM-bioinformatics/ViroConstrictor/compare/v1.5.5...v1.6.0) (2025-12-03)


### Features

* add command line argument for Ampligone's fragment-lookaround-size ([20d2cc5](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/20d2cc58c3943637007c1bea17ea956dcf18f303))
* add pneumoviridae preset ([7d74f6e](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/7d74f6ee8aec5c845d4414e08904f285998ed013))
* added dynamic runtime assigment ([6ec047d](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/6ec047d59882224d3b1437d3ad48f0ee1dca003a))
* added genbank support ([724d3ea](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/724d3ea43936ba482783d43defc5cf02dab695c1))
* added per sample preset disabling in sample sheet with input handling ([56b6cdb](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/56b6cdb965fc83515a5edf9b40b48243a1163912))
* added per-sample setting of AmpliGone's fragment lookaround size in sample sheet with input handling ([9621b73](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/9621b7350c884bfa37ff797f7c73a0afe45f7605))
* added verbose option; add run details to debug messages in general log file ([04878f3](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/04878f39fcd65ccf1d5b536618200e52b92d257e))


### Bug Fixes

* add ambiguity complement bases in clipper.py in order to handle edge cases ([2056c8c](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/2056c8c73486bb01b9853824a25d4b558862f2d5))
* changed default in preset params ([61cb0c8](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/61cb0c868d407ee3ef093c0f2d5d7ff24d7102c2))
* correct the inputs and outputs for all python modules in the match_ref workflow ([d288438](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/d288438fcdc2e6960fbd78e145f4df6002e74a43))
* correctly handle empty primers file in `amplicon_covs.py` when given `--primers NONE` ([2056c8c](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/2056c8c73486bb01b9853824a25d4b558862f2d5))
* correctly output logging of clipper script ([2056c8c](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/2056c8c73486bb01b9853824a25d4b558862f2d5))
* dependency error in consensus.yaml ([61cb0c8](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/61cb0c868d407ee3ef093c0f2d5d7ff24d7102c2))
* do not infer "NA" value in a table as NaN ([325fb8f](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/325fb8f21240f5e47d35293f69e48c5c91a852eb))
* dockerfiles - added adduser command because it is no longer in base ([19f4e28](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/19f4e28e72a353b9bc1edab612de5e4548c5087c))
* enhance genbank parsing with FEATURES column in samplesheet being allowed to inherit from reference genbank file ([54db5a1](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/54db5a18b17ab3cf1b48fd36297981f262eacb34))
* ensure correct feature_type parameter assignment in Translate_AminoAcids rule ([1cd6ae4](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/1cd6ae442cfccba870190e1866fdf743eae2f118))
* fixed amplicon_covs import error ([724d3ea](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/724d3ea43936ba482783d43defc5cf02dab695c1))
* fixed imports in script tests ([6ec047d](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/6ec047d59882224d3b1437d3ad48f0ee1dca003a))
* fixed some copilot PR review spotted errors ([d49b539](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/d49b5399361ae564f81209c4c5a084de3e4f4f91))
* handle empty rows in between rows in sample sheet ([631f633](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/631f633f7a55a245ae548ccdbeb995d76153a103))
* IDSBIO-1334 amplicon covs too strict naming ([d842721](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/d842721bb672ee10f5022f432e880b01f6500c4c))
* improve correction of unidirectional flag logic for illumina platform based on sample input ([9faf0c8](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/9faf0c848ed930bce940f91f53194ad9be9c4aa3))
* **logging:** remove snakemake logrecords from the records_to_emit list when the levelno is lower than the current log level ([a8fe02c](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/a8fe02cbdc2b19645198786b1d05462ef2bcc89c))
* **logging:** set log level for shell and file handlers in ViroConstrictorBaseLogHandler ([a8fe02c](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/a8fe02cbdc2b19645198786b1d05462ef2bcc89c))
* made primer name parsing more stable ([3921e27](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/3921e27af2c4155b4d59285330894b7bbb6cc340))
* match_ref.workflow - regex patterns should be raw strings ([d761cdd](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/d761cddacb252529906231e94a8dfdf962e6250d))
* **presets:** change default feature type to 'gene' and enhance get_preset_parameter function so output AA structures can be found dynamically on a per-preset / per-sample basis ([baede05](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/baede05b90a7369c034e4dbe6675bd525b4e4f68))
* remove own path and remove complexity from amplicon_covs ([e1db2e2](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/e1db2e218dab19f1bd9602116391da8b36b60a62))
* removed "" around snakemake input ([b5ae840](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/b5ae840ff79bf62448c0fd0d33ec39616b7fd9bd))
* replace unnecessary f-strings with normal strings ([d49b539](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/d49b5399361ae564f81209c4c5a084de3e4f4f91))
* resolve path of multiqc config file ([daf0174](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/daf0174dc64ed562ff4b9994c72c4a8cf07db413))
* resolved lambda positional argument error when calling memory requirements function ([63b47ab](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/63b47abdfb1f1dc20476b1e7da5dab99be542261))
* runreport.py - ln=1 was deprecated ([0a87176](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/0a87176fad134a5388612953462bc8ae4e484f5b))
* test_e2e - fixed assertion statements ([d761cdd](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/d761cddacb252529906231e94a8dfdf962e6250d))
* ugly fix for annotated string bug. It edits part of the debugger code to transform snakemake.AnnotatedString types in str types. Could not find an easier way of doing this. ([0398a21](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/0398a21260b4e9cdf6f74dcfacbbf2d8016a1c33))
* update dependencies ([eb70e92](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/eb70e9208e316f5079c20d72658cff486a457c5d))
* update import path for BaseScript in group_refs.py ([881d4fc](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/881d4fc1ae32f52481a6463e23b110f7b4da29a1))
* update script path for filter_bed rule in filter.primers.smk ([2d75a53](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/2d75a53fad6842c134de3390079f4959f946ff45))
* Update ViroConstrictor/workflow/main/scripts/amplicon_covs.py ([d9143e3](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/d9143e392a52e1795b34f8c53e1d5a09ba2c1276))
* update.py - distutils Version classes are deprecated ([0a87176](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/0a87176fad134a5388612953462bc8ae4e484f5b))
* use the new AminoExtract api functions in GFF related scripts ([d288438](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/d288438fcdc2e6960fbd78e145f4df6002e74a43))
* used bool instead of enum for binary state in PrimerInfo class ([cc922cb](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/cc922cb4a22afe15a4868a1133bdf543f057f89c))
* workflow.smk - regex expression should be raw strings ([0a87176](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/0a87176fad134a5388612953462bc8ae4e484f5b))


### Performance Improvements

* slightly reduce size of containers ([a6dfdc9](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/a6dfdc9937fa01695b9869c53ccfb1425be75fe9))


### Dependencies

* also update aminoextract in consensus sequence environment ([bdc7160](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/bdc71605cc57f400c492aad3a083866a3ba555a3))
* change biovalid installation to the bioconda channel and pin version, update pyproject.toml accordingly ([d0fb7c9](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/d0fb7c9cd8de6e7a1f7221a0d9b8d1c6b40a1764))
* fix bug with multiqc where it requires `ps aux` when processing large datasets ([2056c8c](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/2056c8c73486bb01b9853824a25d4b558862f2d5))
* include snakemake executor plugins as base dependencies ([7e3e09d](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/7e3e09df640e4fd99d3ec1b4214f50f39f6cf22e))
* synchronize biopython package through all downstream envs to v1.84 ([ad7f324](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/ad7f324c6e44cbc4172cf7be5cdee04021ab4d18))
* update ampligone to v2.0.2 ([ad7f324](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/ad7f324c6e44cbc4172cf7be5cdee04021ab4d18))
* update bedtools to v2.31.1 ([ad7f324](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/ad7f324c6e44cbc4172cf7be5cdee04021ab4d18))
* update dependencies for all workflow conda environments/containers ([ad7f324](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/ad7f324c6e44cbc4172cf7be5cdee04021ab4d18))
* update dependencies to use AminoExtract version 0.4.0 ([1137c3d](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/1137c3d4d1b33453de846dd4aa09ff545ec75b6a))
* update fastp to v1.0.1 ([ad7f324](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/ad7f324c6e44cbc4172cf7be5cdee04021ab4d18))
* update fastqc to v0.12.1 ([ad7f324](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/ad7f324c6e44cbc4172cf7be5cdee04021ab4d18))
* update minimap2 to v2.28 ([ad7f324](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/ad7f324c6e44cbc4172cf7be5cdee04021ab4d18))
* update multiqc to v1.32 ([ad7f324](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/ad7f324c6e44cbc4172cf7be5cdee04021ab4d18))
* update samtools to v1.21 ([ad7f324](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/ad7f324c6e44cbc4172cf7be5cdee04021ab4d18))
* update TrueConsense to v0.5.2 ([daf0174](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/daf0174dc64ed562ff4b9994c72c4a8cf07db413))
* use AminoExtract version 0.4.1 ([ac75670](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/ac756700f623a3816332c4d53548a17fea0549f0))


### Documentation

* minor changes to presets table and architecture mermaid graphs ([2ad94f5](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/2ad94f5dadfee4dd4df868a2eff6364c606d512c))
* re-wrote most of the docs ([09404d0](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/09404d0cda4446423152b4d0145f30ce3cd488b2))
* re-wrote most of the docs ([42edf4b](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/42edf4b99af61bc06c7b55cd7105fe8fead19bbf))
* render lists as actual lists and minor linguistic changes ([ac546ea](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/ac546ea252ba4d73f7fe5faf7a168d35d9aa101f))

## [1.5.5](https://github.com/RIVM-bioinformatics/ViroConstrictor/compare/v1.5.4...v1.5.5) (2025-07-28)


### Bug Fixes

* **workflow:** fix typo of rulename which leads to unresolvable DAG ([af21030](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/af21030c98737092b6317c1366d63ac380b72fe8))

## [1.5.4](https://github.com/RIVM-bioinformatics/ViroConstrictor/compare/v1.5.3...v1.5.4) (2025-06-26)


### Bug Fixes

* correctly handle situation when --primers NONE is given in combination with --match-ref ([fa64ada](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/fa64ada024ff43fcbe5fe1606d2a32440c3917c8))

## [1.5.3](https://github.com/RIVM-bioinformatics/ViroConstrictor/compare/v1.5.2...v1.5.3) (2025-05-20)


### Features

* add optimized preset for Paramyxoviridae virus family. Combined Measles and Mumps into single preset ([043e6c1](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/043e6c1e50245b5c9427f59c91812c61eb7ed22e))
* add preset for hepatoviridae analysis ([b6ef384](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/b6ef3842f02128bd48611872aa3f81e6520b9edf))


### Bug Fixes

* correction of swapped setting in paramyxoviridae & hepatovirus presets ([de2db8a](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/de2db8a9ed0c168cf32f74f62987a771ecdeb267))


### Documentation

* add "fragmented" mode to documentation ([f6e9bba](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/f6e9bbaf2748b97bfcbddf97db955b29a81c554e))
* add apptainer to the container config doc ([1977437](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/19774375c4c0306f1ada9b8047821ef1700b3af3))
* add authors :) ([1977437](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/19774375c4c0306f1ada9b8047821ef1700b3af3))
* Add the Paramyxoviridae and Hepatovirus presets to the documentation including aliases ([de2db8a](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/de2db8a9ed0c168cf32f74f62987a771ecdeb267))
* spell and grammar checked all docs. ([e431d80](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/e431d8020815dc92fc30aa27819fc4cdd14ddd27))


### Miscellaneous Chores

* empty commit ([8d62c92](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/8d62c92e7ed99d919658265a727d7733ec110ead))

## [1.5.2](https://github.com/RIVM-bioinformatics/ViroConstrictor/compare/v1.5.1...v1.5.2) (2025-04-02)


### Dependencies

* update AmpliGone to version 2.0.1 ([af85890](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/af85890925d667c14fa469800396e47d19a4a781))
* update biopython to version 1.84 ([af85890](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/af85890925d667c14fa469800396e47d19a4a781))
* update fastqc to version 0.12.1 ([af85890](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/af85890925d667c14fa469800396e47d19a4a781))
* update mappy/minimap to version 2.28 ([af85890](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/af85890925d667c14fa469800396e47d19a4a781))
* update openjdk to version 11.0.23 ([af85890](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/af85890925d667c14fa469800396e47d19a4a781))
* update pandas to lenient version 2.2 ([af85890](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/af85890925d667c14fa469800396e47d19a4a781))
* update samtools to version 1.21 ([af85890](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/af85890925d667c14fa469800396e47d19a4a781))

## [1.5.1](https://github.com/RIVM-bioinformatics/ViroConstrictor/compare/v1.5.0...v1.5.1) (2025-03-06)


### Bug Fixes

* clipper filter minimum aligned read length should be a fraction of the reflength instead of a direct input ([3942d99](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/3942d99d607281b85cda29cac6395b57f1d06070))
* suppress snakemake 7 html report licence download error ([5efd010](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/5efd010da8185bd7b2a7b4f9a8050838e739891b))

## [1.5.0](https://github.com/RIVM-bioinformatics/ViroConstrictor/compare/v1.4.6...v1.5.0) (2025-01-24)

### Features

* Add container support to main workflow (f455445)
* Add containerization support to workflow (c29495f)
* Add new global userprofile section for reproducibility method (conda or containers) (e67bb27)


### Bug Fixes

* Properly handle the floating point comparison to check the PRESET_SCORE (f4410bf)
* Create an empty {sample}_primers.bed file when no input primers were given (f7c90e3)
* Do not download containers in dryrun mode (4f6e45d)
* Resolve snakemake AmbiguousRuleException for certain combinations of primer inputs (5cf5fbd)


### Build

* Add base container definitions (836daa3)


### Continuous Integration

* Add GitHub actions workflow for automatic building of containers (41bd864)
* Change upstream registry (ac8a3ac)
* Update GitHub Actions workflows for container publishing (9341f95)
* Add listing of pip modules as a test (b5fed9f)
* Update workflow rules due to deprecation (30be0df)


### Dependencies

* Ensure conda strict channel compatibility (527ad85)
* Limit mamba version to <2.0.0 (941c19b)
* Set flexible pyopenssl dependency to version 24.x.x (0dc3d58)


### Documentation

* Add docstrings to all functions in containers.py (8a5fdfe)
* Enhance documentation with reproducibility settings and multi-reference analysis guidance (71b39ea)

## [1.4.6](https://github.com/RIVM-bioinformatics/ViroConstrictor/compare/v1.4.5...v1.4.6) (2024-10-08)


### Bug Fixes

* properly solve DAG workflow for nonsegmented matched-ref samples ([02a821a](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/02a821a44c3ed3741c65825789ef25ad3e2093c1))

## [1.4.5](https://github.com/RIVM-bioinformatics/ViroConstrictor/compare/v1.4.4...v1.4.5) (2024-09-25)


### Bug Fixes

* incorrect renaming of matched references in non-segmented mode ([c42740e](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/c42740e8a3893f5bd682cbef0edb9a47eefb19fd))

## [1.4.4](https://github.com/RIVM-bioinformatics/ViroConstrictor/compare/v1.4.3...v1.4.4) (2024-09-16)


### Bug Fixes

* remove `intel` channel from selfupdater repositories ([43fe508](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/43fe50812ecb90d693892735431c3fe725159e50))


### Dependencies

* reorder and add details to the conda dependency files for compatibility with strict_channel_priority in conda config ([dc9b993](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/dc9b993fecb3d8840de72cebe82061692ce47269))

## [1.4.3](https://github.com/RIVM-bioinformatics/ViroConstrictor/compare/v1.4.2...v1.4.3) (2024-07-03)


### Dependencies

* Remove 'intel' channel from conda environments ([a2e918a](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/a2e918afc1a39326bb8784c2d79d3bb24731588b))

## [1.4.2](https://github.com/RIVM-bioinformatics/ViroConstrictor/compare/v1.4.1...v1.4.2) (2024-04-03)


### Dependencies

* update AmpliGone to version 1.3.1 ([b6a92c8](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/b6a92c8698b69cc518f73c3d59eda76aea53fbca))

## [1.4.1](https://github.com/RIVM-bioinformatics/ViroConstrictor/compare/v1.4.0...v1.4.1) (2024-03-22)


### Bug Fixes

* resolve pandas FutureWarning for integer casting on a series with a single element ([bb76d06](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/bb76d0650881d3559bc323cfc87fd90ba6746c8b))


### Dependencies

* update AmpliGone to version `1.3.0` ([6f28568](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/6f285680665f3757c6e7b460c24f4cf03ad79627))
* update python and pysam version in Consensus conda environment ([18e4ffb](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/18e4ffbfeadc1da822c8184c740bec21ec08f989))

## [1.4.0](https://github.com/RIVM-bioinformatics/ViroConstrictor/compare/v1.3.1...v1.4.0) (2023-09-18)


### Features

* add "NONE" input for primers and gff's support to the match_ref process ([1581c78](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/1581c7807a1457fdd34b215817fa95a2110c29b5))
* add all the match_ref analysis scripts ([1581c78](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/1581c7807a1457fdd34b215817fa95a2110c29b5))
* add the match_ref snakemake workflow ([1581c78](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/1581c7807a1457fdd34b215817fa95a2110c29b5))
* add the match_reference process ([1581c78](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/1581c7807a1457fdd34b215817fa95a2110c29b5))
* only initiate the match_ref wrapping process whenever it is set on the commandline or through the samplesheet ([1581c78](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/1581c7807a1457fdd34b215817fa95a2110c29b5))


### Bug Fixes

* add fix for unformatted printed message ([98386c1](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/98386c1a70f7163a964110ac6764d3a58ea64811))
* enforce a "None" string instead of a Nonetype used for reference prep ([f4f7490](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/f4f74901ce5584043120335b20b62c02a43362cf))
* ensure "match-ref" and "segmented" columns are always added to samplesheet ([98386c1](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/98386c1a70f7163a964110ac6764d3a58ea64811))
* ensure a string value like "NA" will not be interpreted as `NaN` in pandas.read_csv() ([ecaca5e](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/ecaca5ed908862eb50ef86774919f2dcdebc07b4))
* pad the number in primer names in the amplicon_coverages script so the amplicons are sorted correctly ([cab3a8c](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/cab3a8c289223a666f2ebafc39dc2b9d8c1cab03))
* solve index out of range bug in amplicon_coverages script ([ca95a8a](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/ca95a8a40883ac7ee0787d015a8f980d9f616e30))
* suppress the strict-channel-priorities logmessage from snakemake ([1654891](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/165489190cf215ecb0d9344efa352327d36f94fe))
* workaround to make sure the found segments can properly be 'exploded' in the resulting dataframe ([9db4a42](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/9db4a42b3ab40d4c4d9c4f9b78eeeb81bd3897eb))


### Dependencies

* remove anaconda channel from selfupdating method ([3db1cfb](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/3db1cfb1322135edba53991b2ecfb5216e815e33))
* remove the `anaconda` and `defaults` channels as providers in conda environment recipes, set `nodefaults` ([08ed236](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/08ed236aca43b5aff75beb0381af382542fa30f1))
* replace gffpandas in script to AminoExtract functions ([599da1f](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/599da1f1715d13450698215306a498f3dc67d80c))
* update AminoExtract to version 0.3.1 ([08e4a22](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/08e4a2278a0d359f0dd48da53e030a2a17eab8b6))
* use a forked/patched version of gffpandas instead of the pypi version (circumvent a bug in the original repo) ([08e4a22](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/08e4a2278a0d359f0dd48da53e030a2a17eab8b6))


### Documentation

* update the documentation site to represent the new match-ref and segmented arguments both on commandline and and in samplesheet ([94c5ca0](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/94c5ca07a82b9c12f76b93ab7523adfa998e4585))
* updated docstring for function `check_samplesheet_rows()` ([98386c1](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/98386c1a70f7163a964110ac6764d3a58ea64811))

## [1.3.1](https://github.com/RIVM-bioinformatics/ViroConstrictor/compare/v1.3.0...v1.3.1) (2023-05-26)


### Bug Fixes

* enforce absolute paths for files given through the samplesheet ([40889bf](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/40889bfc3e78160ba96f54c7dc727d662120d5ea))


### Dependencies

* use only snakemake-minimal as workaround for protobuf error in py3.11 ([47c70db](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/47c70db2715bf8933f594040dbff79b74c327e1d))

## [1.3.0](https://github.com/RIVM-bioinformatics/ViroConstrictor/compare/v1.2.6...v1.3.0) (2023-05-03)


### Features

* addition of proper logging functionality for the ViroConstrictor package ([20a952b](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/20a952b0a087c55fda3611e72e4b413cc07a16dc))
* replace snakemake logging output through our own log handler for unified output ([20a952b](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/20a952b0a087c55fda3611e72e4b413cc07a16dc))


### Bug Fixes

* always use the absolute path of given files when parsing from the commandline options ([a5b56e5](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/a5b56e5edbae8dee65324b0d049e427a537c1fcd))
* change certainty values of match_preset_name() function to be explicitly floating points to ensure correct datatypes in downstream check ([a74249c](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/a74249c0b0bd6f3c3adb37c8f45d2225b96f03ba))
* change preset fallback certainty threshold to be slightly more conservative ([7542c46](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/7542c46cc61bd5b66c2a382a1dae122a815b73eb))
* clean-handling of mamba solver mismatching (i.e. CDN mismatch with upstream) ([680730b](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/680730b67e98f1fc1947bbe539e77f6e99630de3))
* ensure variables for preset-warning fallbacks are properly set when no fallbackwarnings have to be logged ([d7ddf0a](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/d7ddf0a97799c260bdbabd28bdf5d79c5deadbf0))
* print the sample key instead of dictionary contents in non-existing reference error log ([f8c7a12](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/f8c7a12707d1311986eaedd4b90187c449f6d3ee))
* Properly show all preset-fallback warnings instead of just the first in the index ([48662db](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/48662dbe86e695b47e0bd2247defb8fc7e365a6c))
* suppress snakemake logging output (workaround for https://github.com/snakemake/snakemake/issues/2089) ([20a952b](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/20a952b0a087c55fda3611e72e4b413cc07a16dc))
* ensure `samples_df` and `samples_dict` always contain the same information ([20a952b](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/20a952b0a087c55fda3611e72e4b413cc07a16dc))


### Documentation

* add citation and DOI information to readme and documentation site ([a2cf117](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/a2cf117b906d947ea750bdf5b1125c3ad18e75ae))
* add explanatory docstrings to various functions ([0e1a368](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/0e1a3686c7e391ec18ca55e0c6818bf8807c1d02))
* fix typo in installation docs ([acc30d9](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/acc30d9d3118498398db3cb54dea52aa08b4ebd4))
* update documentation with dedicated page for presets ([7b2a66f](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/7b2a66f60668ab1c827ba0e796add738fc2293d9))
* updated manual with links to more detailed presets functionality explanation page ([7b2a66f](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/7b2a66f60668ab1c827ba0e796add738fc2293d9))
* updated mkdocs configuration to include new page ([7b2a66f](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/7b2a66f60668ab1c827ba0e796add738fc2293d9))
* updated mkdocs configuration to include the direct copy-button for code blocks ([7b2a66f](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/7b2a66f60668ab1c827ba0e796add738fc2293d9))


### Code Refactoring

* Use a generic (`__main__.py`) top level entry-point instead of the named `ViroConstrictor.py` entrypoint ([20a952b](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/20a952b0a087c55fda3611e72e4b413cc07a16dc))
* re-structure argument parsing functionalities into its own class ([20a952b](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/20a952b0a087c55fda3611e72e4b413cc07a16dc))
* re-structure snakemake run-information and config functionalities into its own class ([20a952b](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/20a952b0a087c55fda3611e72e4b413cc07a16dc))
* remove old shell stdout-coloring method with the rich library ([20a952b](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/20a952b0a087c55fda3611e72e4b413cc07a16dc))
* simplify several functions to ensure a properly defined return ([20a952b](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/20a952b0a087c55fda3611e72e4b413cc07a16dc))
* Use f-strings more consistently for i.e. string concatenation with variables ([20a952b](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/20a952b0a087c55fda3611e72e4b413cc07a16dc))
* add type-hints to all functions ([20a952b](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/20a952b0a087c55fda3611e72e4b413cc07a16dc))


### Miscellaneous Chores

* empty commit ([84a34b9](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/84a34b9ea504bb0a52594803f5f5b9dbfe5665f3))


### Dependencies

* change biopython version to 1.81 ([6bd437a](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/6bd437a17c89c7248e8a8c5ccef03df32b6c4fd5))
* change openpyxl version to 3.1.* ([6bd437a](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/6bd437a17c89c7248e8a8c5ccef03df32b6c4fd5))
* change pandas version to 2.0.* ([6bd437a](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/6bd437a17c89c7248e8a8c5ccef03df32b6c4fd5))
* change snakemake version to 7.25.* in base environment ([6bd437a](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/6bd437a17c89c7248e8a8c5ccef03df32b6c4fd5))
* change urllib3 version to be more lenient (1.26.*) ([6bd437a](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/6bd437a17c89c7248e8a8c5ccef03df32b6c4fd5))
* pin AminoExtract version to 0.2.1 ([6bd437a](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/6bd437a17c89c7248e8a8c5ccef03df32b6c4fd5))
* use more lenient package requirements in setup.py to allow for both py3.10 and py3.11 builds ([944f39d](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/944f39d5e15c2d4c9da45c5ac7da1868d0c85535))

## [1.2.6](https://github.com/RIVM-bioinformatics/ViroConstrictor/compare/v1.2.5...v1.2.6) (2023-03-16)


### Bug Fixes

* add amino-acid names as a tuple to allow for hashing ([9f002c7](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/9f002c7b2814fa15ac47be8f040ee77621adba71))


### Dependencies

* allow for more lenient snakemake versioning &gt;=7.15.2 in base environment ([adc6051](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/adc6051d83e856c89c6e253f984bb9ebec77b83c))

## [1.2.5](https://github.com/RIVM-bioinformatics/ViroConstrictor/compare/v1.2.4...v1.2.5) (2023-03-14)


### Dependencies

* update AmpliGone to version 1.2.1 ([7abaedd](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/7abaedd244db5b19362eb5e6b0ec4ba917cb1a42))
* update minimal mamba version to 1.0.0 ([7ad5fe8](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/7ad5fe852945de7a7abc0be1c80f5a02e032523f))


### Documentation

* add information badges to readme ([53e818f](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/53e818f85eb9b5e08ced64ae3cc6a97b370b822a))
* update installation instructions in docs ([da57e82](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/da57e82ad40d65619d26a91fd4d10d173a9a54a5))

## [1.2.4](https://github.com/RIVM-bioinformatics/ViroConstrictor/compare/v1.2.3...v1.2.4) (2023-03-02)


### Dependencies

* use conda for installing aminoextract instead of pip ([1855470](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/1855470c2ead4510e785ec0f81526259281fe73a))

## [1.2.3](https://github.com/RIVM-bioinformatics/ViroConstrictor/compare/v1.2.2...v1.2.3) (2023-03-01)


### Bug Fixes

* allow input fastq files that contain multiple dots ([07994a2](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/07994a2d5caec1fef067fcef3ba19390bbf3f92a))
* ensure compability between the Influenza preset and illumina short read data ([9a9f1d9](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/9a9f1d985d28064e28cc729d75f8e8039c30b8e9))
* ensure match_preset_name() actually returns all required values ([d49f238](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/d49f2381f430b0f3730ba2991c9360179cac2ff3))
* update permissions for GH-actions workflows ([1f760ba](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/1f760bacf71b069b0e1b07bd043e9d4af0857640))


### Documentation

* add citations file ([1f760ba](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/1f760bacf71b069b0e1b07bd043e9d4af0857640))
* update installation instructions for new environment file and environment creation with Mamba ([84183a7](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/84183a79c4c21543eb69cb0ede0a7923eb54459c))


### Dependencies

* add python-magic version 0.4.27 to environment ([84183a7](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/84183a79c4c21543eb69cb0ede0a7923eb54459c))
* pin version of AminoExtract to 0.2.1 ([84183a7](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/84183a79c4c21543eb69cb0ede0a7923eb54459c))

## [1.2.2](https://github.com/RIVM-bioinformatics/ViroConstrictor/compare/v1.2.1...v1.2.2) (2023-01-24)


### Bug Fixes

* correctly parse presets when using  a samplesheet ([2ca2149](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/2ca21493986cf7480e14f6ada6cdada6cdcf0ba1))
* make the ViroConstrictor main entrypoint less case sensitive ([43b27aa](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/43b27aaa1d8cd95caf01781e964f910649d374fe))

## [1.2.1](https://github.com/RIVM-bioinformatics/ViroConstrictor/compare/v1.2.0...v1.2.1) (2023-01-19)


### Bug Fixes

* correction in argument parser to avoid a NoneType AttributeError when all required information is passed through commandline options ([f2acedf](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/f2acedfb610f87d754fe38303609b528f1afacfc))


### Dependencies

* correct the `rich` version in all base environment files ([f2acedf](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/f2acedfb610f87d754fe38303609b528f1afacfc))

## [1.2.0](https://github.com/RIVM-bioinformatics/ViroConstrictor/compare/v1.1.0...v1.2.0) (2023-01-18)


### Features

* add option to filted for minimum aligned length in `clipper.py` ([d06caf5](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/d06caf5bdb89c39605639df39ee698d54d842a83))
* add option to filter spliced reads in `clipper.py` ([d06caf5](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/d06caf5bdb89c39605639df39ee698d54d842a83))
* add support for working with presets for specific viral targets ([29bfa3c](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/29bfa3c5a1ff0a03ebace707b1c87e298ff1451c))


### Bug Fixes

* add `.gff3` file extension as allowed to features flag ([b4fdf5d](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/b4fdf5dc3cae26639b600adc3a6b7b9a0026240f))
* change Clipper settings to a per-platform config per preset ([3d952bb](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/3d952bb004cea064d67e583a87d72a2a8fa38781))
* change the minimal alignment length in `clipper.py` to be a percentage of the reference length ([68597e8](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/68597e8a5f36e4efc1efaf652563d6e5c4960363))
* correct minimum aligned read length default value ([60b9230](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/60b9230621f0744eaa643aa352d57fed7d96d532))
* Fix extension checking for input files ([6c4f91e](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/6c4f91e812fd1a5c46919413656fa49308dd8ea3))
* set snakemake `keepgoing` to True as a default ([0076d52](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/0076d52c463ba0ad76eb9cdd1d87254e1ffec0c6))
* use the correct columns of input BED file when there are more columns given than necessary ([56b5259](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/56b5259a3fa1803b0b19fb5e9206527bb0b4143c))
* workaround for `Argument list too long` error when parsing complex workflows ([2835c44](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/2835c44de75da7bf59143aeb617ac761ccc8ae56))


### Performance Improvements

* change `make_pickle` rule to localrule ([f41cc01](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/f41cc018d79477b6c9733f6b148ff508b5acb481))
* change threads in local execution mode to improve parallel performance ([869c81a](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/869c81a27c6db56667f56fdaedd3b82f1debaada))


### Dependencies

* add `Rich` as a dependency with version 13.*.* ([e63cac7](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/e63cac778dd4f75c89677eb660bc209fefc3c4da))
* update AmpliGone to version 1.2.0 ([dfb22f1](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/dfb22f10369fa9c4b3b2d89888c70fe101cf0fdf))


### Documentation

* update documentation to include presets functionalities ([adc5265](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/adc526541406363b3cf3afdbf225b71616c9d3a5))

## [1.1.0](https://github.com/RIVM-bioinformatics/ViroConstrictor/compare/v1.0.4...v1.1.0) (2022-12-06)


### Features

* add base of functionality to translate and write amino acid sequences ([25884b7](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/25884b7085d6c4eed616968b749c500494b15e9b))
* translate and extract aminoacids sequences for multiple targets ([657aaff](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/657aaff0928014c7b73826a2392a9b1b49993a80))


### Bug Fixes

* add the missing "FORWARD" and "REVERSE" primer keywords ([74b92c6](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/74b92c648a357059376ba202d7b7d4f0c83ef485))
* change found features to `np.nan` if input reference and features file don't match ([3e812e4](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/3e812e41360e00a176456f1b01895f4084933358))
* circumvent issue where null values cause problems in data translation steps ([6007fca](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/6007fca6f424e8f4499a692d204551b0fc08ca69))
* correctly parse amino acids that have a name with a dot in it ([3e14dd5](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/3e14dd525a75d5a5efe6d59e0b20df6e07cbd53b))
* make sure all required info is present when no gff files are given ([d589989](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/d5899894a9a6b8f2d6cf3e1c0356caa4646bc5d2))


### Documentation

* clarify some sections in the docs ([9bc4171](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/9bc417105773250aa41d6f77ae9517a5545997bd))
* fix broken links ([9bc4171](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/9bc417105773250aa41d6f77ae9517a5545997bd))
* update documentation to match new functionality ([9bc4171](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/9bc417105773250aa41d6f77ae9517a5545997bd))


### Dependencies

* add AminoExtract to dependency list ([acb378e](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/acb378ecb68043d817101d548408636cbfe8adf8))
* change minimal python version to 3.10 ([acb378e](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/acb378ecb68043d817101d548408636cbfe8adf8))
* pin AminoExtract version to v0.2.0 ([cce1cf9](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/cce1cf9fdf71c096dc7d5ca1be5f90be8459c158))
* pin AminoExtract version to v0.2.1 ([c0e869f](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/c0e869f844a1f0c4dfb22cf1ef903a8cf45d10dd))
* update TrueConsense to v0.5.0 and change rules accordingly ([af15b9c](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/af15b9c477540337cf1222c72eaa692e1956798a))
* update TrueConsense version to 0.5.1 ([b1a61e4](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/b1a61e452b219d75948bfbd246e1f7baaa0873a1))

## [1.0.4](https://github.com/RIVM-bioinformatics/ViroConstrictor/compare/v1.0.3...v1.0.4) (2022-11-01)


### Bug Fixes

* cleanly exit vcf_to_tsv.py script when vcf is empty after filtering ([037467b](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/037467b5fa213cffdc45b1f04bd71d83f4a358b9))
* replace bcftools with python script to avoid bcftools shared libraries error ([eefaead](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/eefaeadb07008d00412aafac36579ae3bd30ffaf))


### Dependencies

* remove bcftools as a dependency ([eefaead](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/eefaeadb07008d00412aafac36579ae3bd30ffaf))
* remove the entire Mutations environment as it is no longer necessary ([eefaead](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/eefaeadb07008d00412aafac36579ae3bd30ffaf))

## [1.0.3](https://github.com/RIVM-bioinformatics/ViroConstrictor/compare/v1.0.2...v1.0.3) (2022-10-12)


### Dependencies

* pin snakemake to version 7.15.2 for dependency compatibility ([0914b56](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/0914b561913a0a84c181b9c646d15d34c8409502))
* set snakemake 7.15.2 as a minimal version in setup.py ([0914b56](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/0914b561913a0a84c181b9c646d15d34c8409502))

## [1.0.2](https://github.com/RIVM-bioinformatics/ViroConstrictor/compare/v1.0.1...v1.0.2) (2022-08-15)


### Dependencies

* update AmpliGone to version 1.1.0 ([153719d](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/153719d0aa10862b26024aca4408d1894ab2c552))
* update included snakemake version to 7.12.x ([e59edf9](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/e59edf99de0deeffaf770c3f97883f1b4f1c358a))
* use the bioconda channel to install AmpliGone instead of pip ([0f1fd65](https://github.com/RIVM-bioinformatics/ViroConstrictor/commit/0f1fd65728c1c6ae473f690374a00eda18b6d37a))

### [1.0.1](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/compare/v1.0.0...v1.0.1) (2022-04-26)


### Bug Fixes

* correctly exit amplicon_coverage calc when only one side of an amplicon could be found combined with low coverage samples ([aa77006](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/aa77006992cef2c9f34e7514f3d2c17b73b9526c))
* ensure the split reference sequences are always in uppercase ([fcf034c](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/fcf034ce2bae75731b3e27e69a5e4499a56345b9))

## [1.0.0](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/compare/v0.4.2...v1.0.0) (2022-04-14)


### ⚠ BREAKING CHANGES

* Use a single coverage level during analysis, defaults to 30 ([6d6fc01](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/6d6fc016741cf818817790f66ab85df1ddcd6e63))
* add support for per-sample analysis settings with an input samplesheet in Excel, TSV or CSV format ([a0b726](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/a0b726d544cab9bb61918037e6e0cf04470e0029))
* Add support for BED files in amplicon_covs ([ece3dfa](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/ece3dfa359d1937453a729209a4fbbc8dee1aa44))
* Use BED files as output of AmpliGone ([79fd23](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/79fd2382a60e78000600997b750e703a5c89acc0))
* Update to python version 3.10 in base environments. (minimal required version = 3.8) ([184328](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/18432840a8ec3b2edfebfd101ffe13fb36103cd6))
* Use a single coverage level for TrueConsense ([e5701b9a](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/e5701b9acc50204dec87f5224b9c4193ad57b815))


### Features

* flexible snakemake rule_all based on samplesheet inputs ([93ad0c](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/93ad0cfc55ae46023a1b81129d0f7894aa036703))
* use the snakemake Paramspace for wildcard generation ([94d1738](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/94d1738c1d0a47958f9301261f57338ea544fed8))
* Add remove_adapters_p1 to all ([64c9069](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/64c906939259c43eaaaff3adb2e1c5f8f981e28d))
* Add remove_adapters_p2 and qc_filter to all ([a89560](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/a8956089455320236b328132061f0c5818a77519))
* Add prepare_primers ([6dc7c9](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/6dc7c9321cfbfcd3b92c1f652580573fd2c47e31))
* Add prepare_gff rule ([a0bc0fb](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/a0bc0fb7a00271aaa30c320a62b58f8a0f3e021b))
* Re-add ampligone ([89cc2a0](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/89cc2a0ad9167e1932cf0313a7e309e73971f821))
* Re-add qc_post to multiqc ([45c0d12](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/45c0d124dd993041be49a8433bc0d1f1c6dfbbc0))
* Re-add TrueConsense ([44f36c1](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/44f36c1b922e0da53589c2ac10984f710553d4e3))
* Re-add concat_seq ([ab229de](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/ab229dedb2588f145791da5a72c723301237e8e3))
* Re-add tsv output ([5bf10aa](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/5bf10aa8a14b74f37c7b94c7438f12f2a1bfedf6))
* Re-add breadth of coverage to snakefile ([003b2f](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/003b2f45e5604ed6c8282b0369f243409ba790f1))
* Re-add amplicon_cov ([6f6088](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/6f608826eca77d36abdcc148a0f1d2adea39df34))
* Correctly handle missing inputs in samplesheet ([aa6995](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/aa69950d4e1f3cbec00242ab63ca0ac46226339b))
* re-add illumina support in new workflow structure ([8eb400](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/8eb400bb6712028d6b26c651c251717af8ba60bc))
* add fastP stats to multiqc report ([19627f](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/19627f672473d8214e1ca1502df8dd9cea76271e))


### Bug Fixes

* Add option to not give primers on cli ([e4c054c](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/e4c054c570f670db6ce47bbfea788c4f937cc9ab))
* Check combination of samplesheet and cli args ([e4c054c](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/e4c054c570f670db6ce47bbfea788c4f937cc9ab))
* Fix broken snakefile ([a704c5e](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/a704c5e1ff1bc9fb816d8e5a7f5e56c2663ee3c2))
* Add .bed as filetype for primers ([f16556](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/f165564010da178ff8de86a18792e6f7cc5f5707))
* Make prodigal and prepare_gff unique ([761337](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/7613376599ca116c6b2696898ecf03f635b33124))
* Fix aggregation outputs to rule all ([2269340](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/2269340d95bbb5c821070e22831c56d7118bfff4))
* Fix restructuring of results folder ([76fda46](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/76fda4679c7ee00b258c169573f05a37d1bb6a2d))
* group concatenated files in correct results subfolders ([f3f1d9e](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/f3f1d9e50f65f15b1ea3b5b1f356710eb140f9b7))
* Make sure all post_qc is added in multiqc ([b51f0ce](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/b51f0ceab0460bab1cd5420810dead3699d21b16))
* Remove excess bracket in various log&bench rules ([4a5e96c](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/4a5e96c883fa2e068d5d939a75c9f28aa958510c))
* add logs for calculate_amplicon_cov rule ([fbd2c1c](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/fbd2c1cf5e281375d2c1109443c8cd7d502f862f))
* always use absolute path for input reference, primer and features files ([fa68a7a](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/fa68a7a56c95117915351fa5b101f99c1bdcb3ef))
* always provide absolute path for input fastq files ([0303244](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/030324495165f3e587f5b14da55f633bd7a613fa))
* Raise error if fasta contains invalid char ([5209a7b](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/5209a7b2d385b4779403f8781fe6ab4dd0a3051f))
* Bump minimum snakemake version ([a31b503](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/a31b50332551b4ab767703cc4c0e796e98b752d1))
* change multiqc memory requirements to 'high' ([05b9db3](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/05b9db33e8b26807f9f2c6cf515d7332ca654bd8))
* filter input primer bedfile to only output matching refID bedfile ([c2b9c0c](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/c2b9c0ce3f59d998d6a6dff1520e17cf78fd5efc))


### Code Refactoring

* move version to __init__ file ([e4c054c](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/e4c054c570f670db6ce47bbfea788c4f937cc9ab))
* Remove empty.primers ([bfe0044](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/bfe0044978aa71c871fe74a57cf8f8b7ae1c70fc))
* use the snakemake Paramspace wildcards in the multiQC aggregation ([1204b09](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/1204b09fafb67802caf938353ba9fd79c8e4b1a9))
* add the "extract_refs.py" script to repo ([5a7ea4](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/5a7ea4023a9a8a70ad656c007b9ce90478b2a339))
* move multiqc input function to top of snakefile ([1134bb](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/1134bbf67d6fcb34829bf8b8ad54f5507e8aa3a5))
* use new prepare_primers rule in new structure ([4d06939](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/4d06939eca1c4b76e40ed84ef8906af2bed46d75))
* use new prepare_gffs rule with the new structure ([537d09](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/537d09b596fdf2211a0e02df02b3416e8dbf7305))
* add basic gff filter script ([09ada9](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/09ada979a36765ff2d24b3cd5a2f3af31d8b9838))
* Refactor with bottum up approach ([682cafa](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/682cafa93f74efe137e6de2a0af745dc6e5a226d))
* Rename Sample to sample ([a6c0e48](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/a6c0e4821b48e916c311c8576b330c812fe42784))
* Change to wc_pattern ([214e21](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/214e216c0921f6d069b2b75d4994924d87c8e656))
* Switch wc_pattern to wc_folder ([40d4d2e](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/40d4d2e9aee8edb241cdd646a35170cd9b97dc7b))
* Remove threads: 1 as it is default ([5cd6ca6](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/5cd6ca6b51dab55b6f860648524d3a4b95c6be32))
* Make wildcard notation consistent ([8fe4f7e](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/8fe4f7e3a2eaf057a40d57f44692d14b961455b3))
* Re-ordering of code ([b3d891](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/b3d891be0d53530bfa51a547118061bc3a572c92))
* Remove mincov variable in snakefile ([5bb47a9](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/5bb47a945b33d75cc920c947e59f83eec12aca94))
* change used fonts in the generated report to accomodate fpdf2 ([5bb91d1](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/5bb91d161d28973ac7d61441c78aad066d1d2831))
* merge if statement blocks to singular condition ([7b2882f](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/7b2882fe51a8a87d9e86f72d534745a19a46edc4))
* simplify boolean returns ([615d2cf](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/615d2cf9ad6e84e35d74a6b5bf6437b6cd9461b6))
* use new expression for conditional in samplesheet generator ([8ca585a](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/8ca585a0b4bacb87cc23961457bb5fe3d6b62236))


### Dependencies

* update AmpliGone to version 1.0.1 ([e4c054c](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/e4c054c570f670db6ce47bbfea788c4f937cc9ab))
* update TrueConsense version to v0.4.0 ([e4c054c](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/e4c054c570f670db6ce47bbfea788c4f937cc9ab))
* use Pypi for AmpliGone installation instead of building from source ([e4c054c](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/e4c054c570f670db6ce47bbfea788c4f937cc9ab))
* add python3.7 to the Alignment env for improved env-stability ([d87c746](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/d87c7462af624a01a79acc182db66ec6d0bf719c))
* update snakemake to lenient version 7.3 ([56d8fe](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/56d8fe01f5c872261cd83e8288e2f854af66b579))


### Documentation

* Add information for primer bed support ([e4c054c](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/e4c054c570f670db6ce47bbfea788c4f937cc9ab))
* update docs for multi-target analysis ([e4c054c](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/e4c054c570f670db6ce47bbfea788c4f937cc9ab))
* update documentation for new ViroConstrictor functions ([e4c054c](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/e4c054c570f670db6ce47bbfea788c4f937cc9ab))
* update installation command for python3.10 ([e4c054c](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/e4c054c570f670db6ce47bbfea788c4f937cc9ab))
* update documentation index page with correct required snakemake version ([e4c054c](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/e4c054c570f670db6ce47bbfea788c4f937cc9ab))
* Add explanatory docstring to `CheckSampleProperties` function ([4a552ad](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/4a552ad1ab1e2c3385da904fc760c356c709b1d7))
* Add explanatory docstrings to all functions in `ViroConstrictor.parser` ([b09f4d](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/b09f4dc0d3e7f3db301c871bfd5f924ce91b3b82))
 * Add explanatory docstrings to all functions in `ViroConstrictor.runconfigs` ([2ff8185](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/2ff818572565595fedb7124cdac084fbde8c273d))
* Add explanatory docstrings to all functions in `ViroConstrictor.samplesheet` ([a1e7093](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/a1e7093ff58a1b077131d5b898a7d4e7b85e38e5))
* Add explanatory docstrings to all functions in `ViroConstrictor.userprofile` ([30ec8af8](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/30ec8af8cc2a2493df786e28c35de481460e2cf8))
* Add explanatory docstrings to all functions in `ViroConstrictor.validatefasta` ([c199d7c9](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/c199d7c9b6104cddedce1dff105e9ab713d134cb))
* Add error message for samplesheet duplicates ([efa0761](https://www.github.com/RIVM-bioinformatics/ViroConstrictor/commit/efa0761538a257202a30c0832f3979b790ccb91d))

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
