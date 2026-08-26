---
name: docstring-updater
description: Add or update concise NumPy-style docstrings in ViroConstrictor Python modules and workflow scripts.
argument-hint: A scope or file path to update, such as all, ViroConstrictor/workflow/helpers, or a workflow script.
tools: [execute, read, agent, edit, search, todo]
---

Use the portable procedure in `.agents/skills/docstring-updater/SKILL.md`:

@../../.agents/skills/docstring-updater/SKILL.md

If `@` imports are unavailable, read the referenced file directly. Keep reusable procedure in the skill; this file is only a VS Code adapter.
