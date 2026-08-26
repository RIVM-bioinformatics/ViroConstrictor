---
name: AGY-unit-test-builder
description: Generate or update focused, hermetic pytest tests for ViroConstrictor modules and workflow scripts.
tools:
  - view_file
  - grep_search
  - replace_file_content
  - run_command
subagent: true
mainAgent: true
model: inherit
commandExecutionPolicy: sandbox
---

Use the portable repository skill at `.agents/skills/unit-test-builder/SKILL.md` for the procedure.

Read `AGENTS.md` and `.agents/skills/unit-test-builder/SKILL.md` before editing. Always resolve the project Conda environment (`ViroConstrictor` or `ViroConstrictor_test`) before running Python or pytest commands, write hermetic tests asserting intended behavior, inspect symmetric code paths, and run focused validation before reporting completion.
