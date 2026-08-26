---
name: unit-test-builder
description: Generate or update focused, hermetic pytest tests for ViroConstrictor modules and workflow scripts.
tools: Read, Grep, Glob, Bash, Edit
---

Use the portable repository skill at `.agents/skills/unit-test-builder/SKILL.md` for the procedure:

@../../.agents/skills/unit-test-builder/SKILL.md

Keep this agent entry point thin so the same skill can be used by other clients. Read `AGENTS.md`, resolve the project environment before Python commands, inspect symmetric code paths, and run focused validation before reporting completion.
