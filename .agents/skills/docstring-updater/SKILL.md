---
name: docstring-updater
description: Add or update concise NumPy-style docstrings in ViroConstrictor Python modules and workflow scripts. Use when documenting public APIs without changing runtime behavior.
---

# Docstring Updater

Document public modules, classes, methods, and functions in ViroConstrictor using the repository's NumPy-style conventions. Read `AGENTS.md` and the nearest implementation before editing.

## Procedure

1. Identify the requested scope and public symbols.
2. Read signatures, callers, tests, and existing documentation before describing behavior.
3. Add concise triple-double-quoted docstrings with a summary and applicable `Parameters`, `Returns`, `Raises`, `Notes`, or `Examples` sections.
4. Include parameter and return types when the signature provides them.
5. For `BaseScript` classes, document the input/output paths and CLI arguments at a high level.
6. Change only documentation and harmless formatting; preserve runtime behavior and public signatures.
7. Run a focused syntax, lint, or documentation check after editing and report the result.

Keep `SKILL.md` focused. Put large examples or reference material in separate files only when needed.
