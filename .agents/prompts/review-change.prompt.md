---
name: "review-change"
description: "Review a ViroConstrictor diff, pull request, or review thread for valid behavioral and compatibility risks."
argument-hint: "Provide a PR, diff, issue, review comments, or files to inspect."
---

Use the [deep-dive researcher skill](../skills/deep-dive-researcher/SKILL.md) for evidence gathering.

Review scope: ${input:scope:Provide a PR, diff, issue, review comments, or files to inspect.}

1. Read the relevant implementation, callers, tests, configuration, and history before judging a comment.
2. Trace the complete runtime path, including installed-package and end-user execution where relevant.
3. Classify each finding as valid, invalid, or requiring clarification. Separate confirmed facts from assumptions and distinguish intended design from accidental behavior.
4. Check supported parallel paths, failure behavior, security implications, and the narrowest useful validation command.
5. Apply a Poka-yoke footgun lens to pipeline and configuration changes: look for silent data loss, unsafe defaults, path or quoting mistakes, ordering and wildcard errors, partial or stale outputs, privilege and secret leaks, and mismatches between build, test, publish, and end-user paths. For each material risk, state the failure mode and the smallest guard or validation that would catch it.

Return findings first, ordered by severity, with file links and concrete evidence. Then list open questions, missing tests, and a concise change summary. Do not edit files unless the request explicitly asks for implementation.
