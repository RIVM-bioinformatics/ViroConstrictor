---
name: "test-failure"
description: "Reproduce and fix a ViroConstrictor pytest failure with the correct environment and behavior-focused coverage."
argument-hint: "Paste the pytest failure or name the failing test and module."
---

Use the [unit-test builder skill](../skills/unit-test-builder/SKILL.md) for the testing procedure.

Failure or target: ${input:failure:Paste the pytest failure or name the failing test and module.}

Resolve the project Conda environment before running Python or pytest commands. Read the nearest implementation, fixture, and test before editing. Check expected errors and boundary cases, and inspect symmetric implementation paths and fixtures together. Run the narrowest relevant pytest selection after each substantive edit; widen only after that check passes.

Report the root cause, behavior covered, exact validation command and result, and any dependency or environment blocker. Keep external services, network access, containers, and heavy binaries out of unit tests unless the test explicitly targets that integration boundary.
