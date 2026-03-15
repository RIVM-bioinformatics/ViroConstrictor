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

Typical setup:
```bash
mamba env create -f env.yml -n ViroConstrictor_test
mamba activate ViroConstrictor_test
pip install --upgrade pip
pip install -e .
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
pytest -q tests/unit
```

If failures occur, narrow scope and inspect deeply:
```bash
pytest -q -k "<pattern>" -vv tests/unit
pytest --maxfail=1 --full-trace tests/unit
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
