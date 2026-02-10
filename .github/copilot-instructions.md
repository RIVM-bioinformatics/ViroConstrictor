# ViroConstrictor AI Coding Instructions

## Project Overview
ViroConstrictor is a viral amplicon-sequencing pipeline built on **Snakemake 9.5+**. It processes raw FastQ data (Nanopore, Illumina, IONTorrent) and generates consensus sequences using TrueConsense.

**Key technologies**: Python 3.10+, Snakemake 9.5+, Pandas, Biopython, Pysam, Apptainer/Singularity containers.

## Entry Flow
```
ViroConstrictor/__main__.py 
    → CLIparser (parser.py) 
    → WorkflowConfig (workflow_config.py)
    → WorkflowExecutor (workflow_executor.py) 
    → Snakemake workflows
```

### Key Components
- **[ViroConstrictor/parser.py](ViroConstrictor/parser.py)**: CLI parsing, sample sheet handling, input validation (~1300 lines - the largest module)
- **[ViroConstrictor/workflow_executor.py](ViroConstrictor/workflow_executor.py)**: Snakemake API wrapper using `SnakemakeApi` context manager
- **[ViroConstrictor/workflow_config.py](ViroConstrictor/workflow_config.py)**: Snakemake settings configuration (OutputSettings, ResourceSettings, DeploymentSettings, etc.)
- **[ViroConstrictor/scheduler.py](ViroConstrictor/scheduler.py)**: Multi-scheduler support (LOCAL, SLURM, LSF, DRMAA) with auto-detection
- **[ViroConstrictor/match_ref.py](ViroConstrictor/match_ref.py)**: Match-reference workflow orchestration
- **[ViroConstrictor/genbank.py](ViroConstrictor/genbank.py)**: GenBank file parsing and conversion to FASTA+GFF
- **[ViroConstrictor/samplesheet.py](ViroConstrictor/samplesheet.py)**: Sample detection patterns per platform
- **[ViroConstrictor/runreport.py](ViroConstrictor/runreport.py)**: PDF report generation using FPDF

### Workflow Structure
```
ViroConstrictor/workflow/
├── main/
│   ├── workflow.smk              # Main entrypoint with rule `all`
│   ├── components/
│   │   ├── preparation.references.smk   # Reference FASTA preparation
│   │   ├── preparation.primers.smk      # Primer BED file generation
│   │   ├── preparation.features.smk     # GFF feature preparation + Prodigal fallback
│   │   ├── stats.pre_clean.smk          # FastQC before cleaning
│   │   ├── clean.adapter_removal.smk    # Adapter trimming with Fastp/Clipper
│   │   ├── clean.data_filter.smk        # Quality filtering with Fastp
│   │   ├── clean.primer_removal.smk     # Primer removal with AmpliGone
│   │   ├── stats.post_clean.smk         # FastQC after cleaning + MultiQC
│   │   ├── results.sequences.smk        # Alignment, consensus, amino acid extraction
│   │   ├── results.reporting_metrics.smk # Coverage, mutations, breadth of coverage
│   │   ├── results.concatenations.smk   # Combine per-sample results
│   │   └── results.combined.smk         # Aggregate results across samples
│   ├── configs/
│   │   └── multiqc.yaml          # MultiQC configuration
│   └── scripts/                  # Python scripts called by rules
│       ├── prepare_refs.py       # Reference sequence preparation
│       ├── extract_gff.py        # GFF feature extraction
│       ├── filter_bed_input.py   # BED file filtering by reference
│       ├── clipper.py            # BAM read filtering/clipping
│       ├── amplicon_covs.py      # Amplicon coverage calculation
│       ├── concat_amplicon_covs.py # Concatenate coverage files
│       ├── boc.py                # Breadth of coverage calculation
│       ├── vcf_to_tsv.py         # VCF to TSV conversion
│       ├── group_aminoacids.py   # Group amino acid sequences
│       ├── combine_fasta.py      # Combine FASTA files
│       ├── combine_tabular.py    # Combine TSV/CSV files
│       ├── aggregate_combined_files.py # Aggregate combined results
│       └── fastqc.sh             # FastQC wrapper script
├── match_ref/
│   ├── workflow.smk              # Match-reference entrypoint
│   ├── components/
│   │   ├── preparation.references.smk   # Filter references by segment
│   │   ├── selection.alignment.smk      # Minimap2 alignment for selection
│   │   ├── selection.reference.smk      # Best reference selection + grouping
│   │   ├── filter.features.smk          # GFF filtering by selected reference
│   │   └── filter.primers.smk           # Primer BED filtering by reference
│   └── scripts/
│       ├── filter_references.py         # Filter FASTA by wildcard segment
│       ├── count_mapped_reads.py        # Count reads per reference
│       ├── filter_best_matching_ref.py  # Select best matching reference
│       ├── group_refs.py                # Group references and stats
│       ├── filter_gff.py                # Filter GFF by reference data
│       └── filter_bed.py                # Filter BED by reference ID
├── helpers/
│   ├── base_script_class.py      # BaseScript ABC for all workflow scripts
│   ├── containers.py             # Container hash calculation and download
│   ├── directories.py            # Path constants for all output directories
│   ├── generic_workflow_methods.py # Shared workflow helper functions
│   ├── presets.py                # Preset matching and parameter retrieval
│   ├── preset_params.json        # Preset configurations (SARSCOV2, INFLUENZA, etc.)
│   └── preset_aliases.json       # Alias mappings for fuzzy preset matching
└── envs/                         # Conda environment definitions
    ├── Alignment.yaml            # Minimap2, Samtools
    ├── Clean.yaml                # Fastp, AmpliGone
    ├── Consensus.yaml            # TrueConsense
    ├── core_scripts.yaml         # Pandas, Biopython, Pysam
    ├── mr_scripts.yaml           # Match-ref specific dependencies
    └── ORF_analysis.yaml         # Prodigal, AminoExtract
```

## Workflow Script Pattern

All workflow scripts inherit from `BaseScript` in [helpers/base_script_class.py](ViroConstrictor/workflow/helpers/base_script_class.py):

```python
from helpers.base_script_class import BaseScript  # type: ignore[import]

class MyScript(BaseScript):
    def __init__(self, input: Path | str, output: Path | str, custom_arg: str) -> None:
        super().__init__(input, output)
        self.custom_arg = custom_arg

    @classmethod
    def add_arguments(cls, parser: ArgumentParser) -> None:
        super().add_arguments(parser)
        parser.add_argument("--custom_arg", required=True)

    def run(self) -> None:
        # Implementation here
        pass

if __name__ == "__main__":
    MyScript.main()
```

### 2. Script Invocation Pattern in Snakemake Rules
```snakemake
params:
    script="-m main.scripts.my_script",
    pythonpath=f'{Path(workflow.basedir).parent}'
shell:
    """
    PYTHONPATH={params.pythonpath} \
    python {params.script} \
    --input {input} \
    --output {output} \
    --custom_arg {params.custom_arg}
    """
```

## Conventions

### Preset System
Presets auto-configure tool parameters (Minimap2, Fastp, AmpliGone) per virus type. Use fuzzy matching in `match_preset_name()` with ≥40% similarity threshold. Always fall back to "DEFAULT" preset.

Preset parameters are retrieved via:
```python
from ViroConstrictor.workflow.helpers.presets import get_preset_parameter

params:
    alignmentpreset=lambda wc: get_preset_parameter(
        preset_name=SAMPLES[wc.sample]["PRESET"],
        parameter_name=f"AmpliGone_AlignmentPreset_{config['platform']}",
    ),
```

### Memory/Runtime Functions (workflow.smk)
```python
def low_memory_job(wildcards, threads, attempt):  # attempt * threads * 1 * 1000
def medium_memory_job(...)                        # attempt * threads * 2 * 1000
def high_memory_job(...)                          # attempt * threads * 4 * 1000
# Runtime scales quadratically: attempt * attempt * [2|15|30] minutes

# Local execution caps memory at config["max_local_mem"]
```
These functions are duplicated in both `main/workflow.smk` and `match_ref/workflow.smk` and in tests - keep all in sync when modifying.

### Directory Constants
Import from `ViroConstrictor.workflow.helpers.directories`:
```python
# Key directories
res = "results/"
datadir = "data/"
cln = "cleaned_fastq/"
aln = "alignment/"
prim = "primers/"
features = "features/"
matchref = "match_ref_process/"
logdir = "logs/"

# Subdirectories
qc_pre = "FastQC_pretrim/"
qc_post = "FastQC_posttrim/"
amino = "aminoacids/"
cons = "consensus/"
covs = "coverages/"
```

### Wildcard Pattern (Paramspace)
Both workflows use Snakemake's `Paramspace` for wildcard management:
```python
# Main workflow
p_space = Paramspace(samples_df[["Virus", "RefID", "sample"]], filename_params=["sample"])
wc_folder = "/".join(p_space.wildcard_pattern.split("/")[:-1]) + "/"
# Results in: "Virus~{Virus}/RefID~{RefID}/"

# Match-ref workflow
p_space = Paramspace(samples_df[["Virus", "segment", "sample"]], filename_params=["sample"])
```

### Container System
Containers are identified by hash of their environment YAML:
```python
from ViroConstrictor.workflow.helpers.containers import get_hash

container:
    f"{container_base_path}/viroconstrictor_clean_{get_hash('Clean')}.sif"
```

Container base path is set via `workflow.deployment_settings.apptainer_prefix`.
Usage of containers requires Apptainer/Singularity to be installed on the system. Docker is only used for building containers and cannot be used for execution.

### Logging Architecture
Use `from ViroConstrictor.logging import log` in core python modules of ViroConstrictor - supports Rich markup:
```python
log.info("[green]Success[/green]")
log.warning("[yellow]Warning message[/yellow]")
log.error("[bold red]Error occurred[/bold red]")
```

ViroConstrictor has a unified logging system that captures both Python application logs and Snakemake workflow logs:

```
snakemake_logger_plugin_viroconstrictor/
├── __init__.py    # LogHandler class (Snakemake plugin interface)
└── handler.py     # ViroConstrictorLogHandler (bridges to base handler)

ViroConstrictor/logging.py
└── ViroConstrictorBaseLogHandler  # Core handler with Rich formatting + file output
```

**How it works**:
1. `ViroConstrictorBaseLogHandler` in `logging.py` handles all log output (console + file)
2. `snakemake_logger_plugin_viroconstrictor` is a Snakemake logger plugin that reroutes Snakemake's internal logs
3. `ViroConstrictorLogHandler` in the plugin extends the base handler as a Snakemake plugin hook
4. Both Python app and Snakemake logs emit to the same handlers with consistent formatting

## Two-Stage Workflow Execution

ViroConstrictor can run in two modes controlled by `vc_stage`:

1. **Match-Reference Stage (`MR`)**: When `MATCH-REF` column is True in samples
   - Aligns reads to multi-reference FASTA
   - Selects best matching reference per sample
   - Outputs: `match_ref_results.pkl` with selected references

2. **Main Stage (`MAIN`)**: Standard analysis
   - Uses single reference (or match-ref selected reference)
   - Full cleaning → alignment → consensus pipeline

The stage is set in `WorkflowConfig`:
```python
if self.vc_stage == "MR":
    self.samplesheetfilename = "samples_mr"
    self.workflow_file = parsed_inputs.match_ref_snakefile
if self.vc_stage == "MAIN":
    self.samplesheetfilename = "samples_main"
    self.workflow_file = parsed_inputs.snakefile
```

## Development Commands

### Testing
```bash
pytest tests/unit/           # Fast unit tests
pytest tests/e2e/            # End-to-end (downloads test data from ENA)
pytest -k "test_scheduler"   # Run specific test
pytest -k "test_amplicon"    # Run amplicon-related tests
```

Test data is located in:
- `tests/unit/data/` - Static test files
- `tests/e2e/data/` - Reference files for E2E tests

### Code Style
- Line length: 150 (black, isort, pylint)
- Python 3.10+ required
- Type hints expected on all new functions
- Docstrings: NumPy style with Parameters/Returns sections

### Container Development
```bash
# Build all containers locally (requires Docker)
python build_local_containers.py

# Build containers that do not exist yet in upstream registry (runs in github actions)
python containers/build_containers.py

# Pull published containers
python containers/pull_published_containers.py
```

Dockerfiles follow naming: `{environment_name}.dockerfile` (e.g., `Clean.dockerfile` for `envs/Clean.yaml`)

## File Patterns

### GenBank Support
`GenBank` class in [genbank.py](ViroConstrictor/genbank.py) auto-splits `.gb`/`.gbk` files into `.fasta` + `.gff` for the workflow. Check with `GenBank.is_genbank(path)`.

### Sample Detection
Platform-specific regex patterns in [samplesheet.py](ViroConstrictor/samplesheet.py):
- **Illumina**: `(.*)(_|\.)R?(1|2)(?:_.*\.|\..*\.|\.)f(ast)?q(\.gz)?`
- **Nanopore/IONTorrent**: `(.*)\.f(ast)?q(\.gz)?`

### Primer File Formats
Supported formats (see [docs/input-formatting.md](docs/input-formatting.md)):
- **FASTA**: Primer sequences with naming `{name}_{number}_{LEFT|RIGHT}`
- **BED**: 6-column format with reference ID matching FASTA header

### Adding New Workflow Rules
1. Create component file in `workflow/main/components/`
2. Include in `workflow.smk`: `include: workflow.source_path("components/yourfile.smk")`
3. Add outputs to `construct_all_rule()` function
4. If new script needed:
   - Create in `workflow/main/scripts/` inheriting `BaseScript`
   - Add unit test in `tests/unit/scripts/`

### Adding New Workflow Scripts
1. Inherit from `BaseScript`
2. Implement `__init__`, `add_arguments`, and `run` methods
3. Use `self.input` and `self.output` from base class
4. Call via `python -m main.scripts.scriptname` pattern
5. Add corresponding unit test

## External Dependencies
- **Snakemake executor plugins**: `slurm`, `lsf`, `drmaa`
- **Bioinformatics tools**: AminoExtract, Minimap2, Samtools, Fastp, AmpliGone, TrueConsense, Prodigal
- **Python packages**: pandas, biopython, pysam, rich, fpdf
- **Validation**: biovalid, bcbio-gff
- **Containers**: Apptainer/Singularity preferred, Docker for building

## Configuration Files
- `~/.ViroConstrictor_defaultprofile.ini` - User configuration (computing mode, scheduler)
- `{workdir}/config.yaml` - Runtime Snakemake config (generated)
- `{workdir}/samples_main.yaml` / `samples_mr.yaml` - Sample sheets (generated)

## Snakemake API Integration
ViroConstrictor uses Snakemake 9.5+ Python API via `SnakemakeApi` context manager:
```python
from snakemake.api import SnakemakeApi

with SnakemakeApi(settings) as snakemake_api:
    workflow_api = snakemake_api.workflow(...)
    dag_api = workflow_api.dag(...)
    dag_api.execute_workflow(...)
```

Key settings classes from `snakemake.api`:
- `ConfigSettings`, `DAGSettings`, `DeploymentSettings`
- `ExecutionSettings`, `OutputSettings`, `ResourceSettings`
- `SchedulingSettings`, `StorageSettings`, `WorkflowSettings`
