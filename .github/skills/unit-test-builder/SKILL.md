---
name: unit-test-builder
description: Generate or update pytest unit tests for ViroConstrictor modules. Use when asked to add test coverage for specific files or the whole codebase, including happy paths, error paths, and edge cases suitable for CI.
compatibility: Requires local Python test environment and repository access; should avoid network and heavy external binaries in unit tests.
metadata:
  source: migrated-from-.github/agents/UnitTestBuilder.agent.md
---

## Purpose
Create and maintain fast, deterministic pytest unit tests for ViroConstrictor.

## Inputs
- A target module path, list of files, or `all`.
- Optional coverage expectations (for example, `>=80%` per touched file).

## Environment Setup
Prefer repository guidance:
- `docs/installation.md`
- `.github/copilot-instructions.md`

Environment selection is mandatory and must be enforced on every command.

Rules:
- Never run tests, package installs, or Python commands in `base`.
- Always use the project environment for this repository.
- Prefer `conda run -n <env> ...` in automated flows to avoid activation drift.

Environment resolution order:
1. Use `ViroConstrictor` if it exists.
2. Else create and use `ViroConstrictor_test`.

Preflight checks:
```bash
conda env list
conda run -n ViroConstrictor python -V
```

If `ViroConstrictor` does not exist:
```bash
mamba env create -f env.yml -n ViroConstrictor_test
conda run -n ViroConstrictor_test python -V
```

Typical setup:
```bash
mamba env create -f env.yml -n ViroConstrictor_test
mamba activate ViroConstrictor_test
pip install --upgrade pip
pip install -e .
```

Recommended non-interactive setup (preferred for agents):
```bash
# Resolve environment once and reuse it for all commands
if conda env list | grep -qE '^ViroConstrictor\s'; then
  VC_ENV=ViroConstrictor
else
  VC_ENV=ViroConstrictor_test
fi

conda run -n "$VC_ENV" pip install --upgrade pip
conda run -n "$VC_ENV" pip install -e .
```

## Test Generation Rules
1. Use pytest idioms (`tmp_path`, `monkeypatch`, `caplog`, parametrization).
2. Name files `tests/unit/test_<module>.py`.
3. Keep tests hermetic: no network, no real container runtime, no heavy binaries.
4. Mock external tools (`apptainer`, `samtools`, `minimap2`, Snakemake subprocess paths).
5. Cover:
   - Happy path behavior.
   - Error handling and invalid input.
   - Edge cases (empty input, missing columns, malformed formats).
6. For workflow scripts, prefer class instantiation and `run()` over subprocess invocation when practical.
7. Keep fixtures concise and reusable.
8. Do not assume the implemented code is correct, test for intended behavior instead of replicating implementation details. It's better to have a test that fails due to an incorrect implementation than a test that passes because it assumes the implementation is correct.

## Module-Specific Expectations
- Parsers (`parser.py`, `samplesheet.py`): success cases, invalid formats, missing columns, useful errors.
- Workflow scripts (`workflow/.../scripts/*.py`): output creation and minimal expected content.
- Helpers: boundaries, parameter retrieval, deterministic utility behavior.
- `genbank.py`: valid GenBank conversion and non-GenBank handling.
- `presets.py`: fuzzy matching, `DEFAULT` fallback, parameter retrieval.
- `containers.py`: hash determinism and mocked download flows.
- `logging.py`: caplog-based behavior and edge-message safety.

## Execution and Failure Handling
Run unit tests after generation:
```bash
conda run -n "$VC_ENV" pytest -q tests/unit
```

If failures occur, narrow scope and inspect deeply:
```bash
conda run -n "$VC_ENV" pytest -q -k "<pattern>" -vv tests/unit
conda run -n "$VC_ENV" pytest --maxfail=1 --full-trace tests/unit
```

Report failures with:
- Failing file and test name.
- Exception type and key traceback line(s).
- Likely root cause.
- Suggested remediation (mocking, fixture/input fix, version mismatch correction).

If a test inherently requires integration tooling, either:
- Refactor to proper unit isolation with mocks, or
- Mark/skip it explicitly as integration-only (only when unavoidable).

## Output
- New/updated files under `tests/unit/`.
- Short coverage summary of branches/scenarios added.
- Notes on any intentionally skipped integration-dependent behavior.
