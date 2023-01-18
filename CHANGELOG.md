# Changelog

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
