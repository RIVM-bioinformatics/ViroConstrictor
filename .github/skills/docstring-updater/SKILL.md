---
name: docstring-updater
description: Add or update NumPy-style docstrings in ViroConstrictor Python modules and workflow scripts. Use when asked to document missing/incomplete docstrings, standardize docstring style, or improve API documentation quality without changing runtime behavior.
compatibility: Designed for repository-local edits in VS Code/Copilot-style coding agents.
metadata:
  source: migrated-from-.github/agents/DocstringUpdater.agent.md
---

## Purpose
Add or update docstrings so public modules, classes, methods, and functions follow NumPy docstring conventions used in this repository.

## Inputs
- A scope string such as `all`, `ViroConstrictor/workflow/helpers`, or `workflow/main/scripts/prepare_refs.py`.
- Optional constraints from the user (for example, "only workflow scripts" or "no signature changes").

## Scope Rules
- Primary targets:
  - `ViroConstrictor/`
  - `ViroConstrictor/workflow/`
  - Main entrypoints and scripts invoked by Snakemake
- Exclude auto-generated files, vendored third-party code, and tests unless explicitly requested.
- Keep behavior unchanged. Limit edits to docstrings and harmless formatting.

## Docstring Standard
Use NumPy style with these requirements:
- Triple double quotes (`"""`).
- One-line summary, then optional detail paragraph.
- Include sections when applicable: `Parameters`, `Returns`, `Raises`, `Notes`, `Examples`.
- Include parameter and return types in docstring sections.
- Respect repository line length guidance (<= 150 chars).
- Keep language concise, technical, and neutral.

For `BaseScript`-derived workflow scripts, document:
- Class purpose.
- Expected input/output path formats.
- CLI argument intent at a high level.

## Quality Checklist
- Every public symbol has a meaningful summary.
- `Parameters` entries use `name : type` and short descriptions.
- `Returns` is present for non-`None` returns.
- `Raises` documents intentional exceptions.
- `Examples` only for non-trivial behavior.
- Preserve or improve type hints in signatures when needed, but do not refactor logic.

## Validation
After edits:
1. Run project checks if available (`black`, `isort`, lint checks) for touched files.
2. Ensure no functional code paths were altered.
3. Keep diff focused on documentation and formatting.

## Output
- Updated Python files with NumPy-style docstrings.
- Brief summary of what was documented and any intentionally skipped files.
- Suggested commit message format: `docs(<scope>): add NumPy docstrings to <target>`.
