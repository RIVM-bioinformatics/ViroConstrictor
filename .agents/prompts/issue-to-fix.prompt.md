---
name: "issue-to-fix"
description: "Investigate a ViroConstrictor issue, failure, or review comment and implement a behavior-backed fix."
argument-hint: "Describe the issue, failure, review comment, or affected files."
---

Use the [deep-dive researcher skill](../skills/deep-dive-researcher/SKILL.md) for the investigation. Use the [unit-test builder skill](../skills/unit-test-builder/SKILL.md) only when the issue requires creating or changing unit tests; an implementation fix does not automatically require test work.

Request: ${input:request:Describe the issue, failure, review comment, or affected files.}

1. Write a compact reproduction brief: trigger, expected behavior, observed behavior, relevant environment, smallest reproducer (or why unavailable), and current evidence.
2. Identify the nearest code path that controls the behavior and read its closest tests or callers.
3. State the intended behavior, one falsifiable hypothesis, and one focused check before editing.
4. Implement the smallest coherent fix. Add or update a hermetic behavior test only when the fix introduces or changes a meaningful testable contract, exposes a regression worth pinning down, or the request explicitly asks for test coverage.
5. Run the narrowest relevant validation immediately, repair local failures, and only then widen the check.

Report a lightweight evidence receipt: diagnosis, trigger and expected/observed behavior, changed files, validation command and result, remaining risks, and any environmental limitation. Distinguish observed results from assumptions, omit unavailable fields rather than inventing them, and do not claim completion when validation was not run.
