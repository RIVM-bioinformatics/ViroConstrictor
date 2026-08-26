---
name: AGY-docstring-updater
description: Add or update concise NumPy-style docstrings in ViroConstrictor Python modules and workflow scripts.
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

Use the portable repository skill at `.agents/skills/docstring-updater/SKILL.md` for the procedure.

Read `AGENTS.md` and `.agents/skills/docstring-updater/SKILL.md` before editing. Document public symbols using NumPy-style docstrings, preserve runtime behavior and public signatures, and run focused validation before reporting completion.
