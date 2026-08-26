---
name: "Python and Test Workflow"
description: "Use for Python implementation and pytest changes in ViroConstrictor; resolve the project environment and validate intended behavior."
applyTo: "**/*.py"
---

- Before running Python, pytest, package, or build commands, use the existing `ViroConstrictor` Conda environment; otherwise use `ViroConstrictor_test` from `env.yml`. Verify the interpreter and never silently run project commands in `base`.
- For a test failure, read the nearest implementation, fixture, and test before editing. State the intended behavior and one focused check that could disprove the diagnosis.
- Run the narrowest relevant pytest selection immediately after an edit. Widen only after the focused check passes or after repairing a failure in the same slice.
- Test behavior rather than implementation details, including expected errors, malformed input, empty input, and important boundaries.
- Inspect symmetric branches and fixtures together, especially LOCAL, SLURM, LSF, and DRMAA scheduler paths. Keep tests hermetic with `tmp_path`, mocks, and stubs for optional heavy dependencies.
