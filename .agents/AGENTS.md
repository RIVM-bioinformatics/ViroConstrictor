# ViroConstrictor Agent Guide

This is the canonical, client-neutral guide for coding agents working in ViroConstrictor. It is intentionally plain Markdown so it can be used by Codex, Claude Code, Gemini, Antigravity, VS Code, and other agents. Explicit user instructions take precedence.

## Project

ViroConstrictor is a viral amplicon-sequencing pipeline built with Python 3.10+, Snakemake 9.5+, Pandas, Biopython, Pysam, and Apptainer/Singularity containers. It processes Nanopore, Illumina, and IONTorrent FastQ data and generates consensus sequences with TrueConsense.

The main execution flow is:

```
ViroConstrictor/__main__.py -> parser.py -> workflow_config.py -> workflow_executor.py -> Snakemake workflows
```

Important areas:

- `ViroConstrictor/parser.py`: CLI parsing, sample-sheet handling, and input validation.
- `ViroConstrictor/workflow_config.py`: Snakemake settings and stage configuration.
- `ViroConstrictor/workflow_executor.py`: Snakemake 9.5+ API integration.
- `ViroConstrictor/scheduler.py`: LOCAL, SLURM, LSF, and DRMAA scheduler selection.
- `ViroConstrictor/genbank.py`: GenBank parsing and FASTA/GFF conversion.
- `ViroConstrictor/samplesheet.py`: platform-specific sample detection.
- `ViroConstrictor/workflow/`: main and match-reference Snakemake workflows, scripts, helpers, and Conda environments.

## Environment And Commands

Before running Python, pytest, pip, package-installation, or project build commands:

1. Detect the available environment manager (`conda` or `mamba`) and inspect environments.
2. Use the existing `ViroConstrictor` environment when available.
3. Otherwise use or create `ViroConstrictor_test` from `env.yml`.
4. Verify the selected interpreter before executing project commands.
5. Prefer non-interactive environment-prefixed commands such as `conda run -n "$VC_ENV" ...` so commands do not silently run in `base` or another interpreter.
6. If no suitable environment manager or environment is available, report the blocker instead of running project commands in `base`.

A POSIX-shell example is:

```bash
if conda env list | grep -qE '^ViroConstrictor\s'; then
  VC_ENV=ViroConstrictor
else
  VC_ENV=ViroConstrictor_test
fi
conda run -n "$VC_ENV" python -V
```

Translate the equivalent steps for the host shell when necessary. Apptainer/Singularity is preferred for workflow execution; Docker is used for building images, not for running the pipeline.

## Testing And Validation

Use focused checks first, then widen only when needed:

```bash
conda run -n "$VC_ENV" pytest -q tests/unit/test_<target>.py
conda run -n "$VC_ENV" pytest -q tests/unit
conda run -n "$VC_ENV" pytest -q tests/e2e
```

For every test change or test-failure fix:

- Test the intended behavior, not merely the current implementation.
- Run the narrowest relevant test after editing and before further unrelated investigation.
- Inspect all parallel or symmetric branches, fixtures, parametrizations, and platform cases. For example, scheduler tests must cover LOCAL, SLURM, LSF, and DRMAA paths where the implementation supports them.
- Check both success and failure behavior, including missing inputs, malformed formats, and unavailable optional dependencies.
- Do not report a fix as complete until the focused check passes or the remaining failure is clearly reported with its root cause.

Keep unit tests hermetic: avoid network access, real container runtimes, and heavy bioinformatics binaries. Mock external tools and use `tmp_path` for temporary files.

## Python And Workflow Conventions

- Add type hints to new Python functions and use NumPy-style docstrings for public APIs.
- Keep lines at or below 150 characters.
- Import logging from `ViroConstrictor.logging` in core modules.
- Use directory constants from `ViroConstrictor.workflow.helpers.directories`.
- Workflow scripts should inherit from `BaseScript`, implement `__init__`, `add_arguments`, and `run`, and be invoked with the `python -m main.scripts.<name>` pattern.
- Add a unit test under `tests/unit/` for new workflow scripts and meaningful behavior changes when the fix introduces or changes a testable contract, exposes a regression worth pinning down, or the task requests coverage. Do not add or modify unit tests automatically for every issue fix.
- Preserve the duplicated memory/runtime helpers in both workflow entrypoints and their tests when changing those functions.
- Use the existing preset system and retain the `DEFAULT` fallback for unknown presets.

The workflow has two stages: match-reference (`MR`) selects a reference from multiple candidates, and main (`MAIN`) performs cleaning, alignment, consensus, reporting, and aggregation.

## Change Workflow

1. Identify the smallest code path that controls the requested behavior.
2. Read the nearest implementation and relevant test or call site.
3. State a falsifiable local hypothesis and choose a cheap check that could disprove it.
4. Make the smallest coherent edit.
5. Run the focused validation immediately.
6. Repair local failures and rerun the same check before widening scope.
7. Summarize changed files, validation commands, and any unresolved environmental limitation.

Do not revert unrelated work in a dirty worktree. Do not commit or create branches unless explicitly requested.

## Agent Customization

This file is the entry point for repository work. Details about organizing skills, prompts, instructions, custom agents, and client adapters live in `.agents/CUSTOMIZATION.md`.
