---
name: unit-test-builder
description: Generate or update focused, hermetic pytest tests for Python modules and workflow scripts. Use when adding coverage, fixing test failures, or validating behavior and edge cases.
---

# Unit Test Builder

Create fast, deterministic pytest tests for ViroConstrictor. Read `AGENTS.md` and the nearest implementation before editing.

## Environment

Before any Python, pytest, pip, or package command, resolve the project environment:

1. Prefer the existing `ViroConstrictor` Conda environment.
2. Otherwise use or create `ViroConstrictor_test` from `env.yml`.
3. Verify the selected interpreter with `python -V`.
4. Prefix commands with the selected environment, preferably `conda run -n "$VC_ENV" ...`. Never run project tests or installs in `base`.

If the required environment cannot be selected, report the blocker instead of silently using another interpreter.

## Procedure

1. Identify the public behavior and the smallest relevant test file.
2. Add or update tests for success, expected errors, malformed or empty input, and important boundaries.
3. For every parallel implementation path, inspect and test the matching fixtures and parametrizations. Do not stop after fixing only the first reported error.
4. Keep tests hermetic: use `tmp_path`, mock network and external executables, and stub optional heavy imports when needed.
5. Run the narrowest relevant test immediately after editing.
6. If it fails, repair the same slice and rerun it before widening the test scope.
7. Report the command, result, and any environmental limitation.

## Project Cases

- Parsers and sample sheets: valid input, missing columns, invalid formats, and useful errors.
- Workflow scripts: instantiate `BaseScript` classes and call `run()` where practical; verify output files and minimal content.
- GenBank: valid conversion, plain FASTA handling, malformed input, and output paths.
- Scheduler: LOCAL, SLURM, LSF, and DRMAA branches supported by the implementation, including symmetric fixtures.
- Helpers and presets: boundary values, deterministic behavior, parameter retrieval, and `DEFAULT` fallback.
- Containers: deterministic hashes and mocked download/build calls.

Tests should assert intended behavior rather than mirror implementation details.
