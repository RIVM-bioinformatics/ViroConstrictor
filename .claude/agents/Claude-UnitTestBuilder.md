---
name: Claude-UnitTestBuilder
description: Generate or update focused, hermetic pytest tests for ViroConstrictor modules and workflow scripts.
tools: Read, Grep, Glob, Bash, Edit
---

Use the portable repository skill at `.agents/skills/unit-test-builder/SKILL.md` for the procedure.

Read `AGENTS.md` and `.agents/skills/unit-test-builder/SKILL.md` before editing. Always resolve the project Conda environment (`ViroConstrictor` or `ViroConstrictor_test`) before running Python or pytest commands, write hermetic tests asserting intended behavior, inspect symmetric code paths, and run focused validation before reporting completion.
