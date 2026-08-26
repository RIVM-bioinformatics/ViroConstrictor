---
name: "issue-to-fix"
description: "Investigate a ViroConstrictor issue, failure, or review comment and implement a behavior-backed fix."
argument-hint: "Describe the issue, failure, review comment, or affected files."
---

Use the [deep-dive researcher skill](../skills/deep-dive-researcher/SKILL.md) for the investigation. Use the [unit-test builder skill](../skills/unit-test-builder/SKILL.md) only when the issue requires creating or changing unit tests; an implementation fix does not automatically require test work.

Request: ${input:request:Describe the issue, failure, review comment, or affected files.}

1. Identify the nearest code path that controls the behavior and read its closest tests or callers.
2. State the intended behavior, one falsifiable hypothesis, and one focused check before editing.
3. Implement the smallest coherent fix. Add or update a hermetic behavior test only when the fix introduces or changes a meaningful testable contract, exposes a regression worth pinning down, or the request explicitly asks for test coverage.
4. Run the narrowest relevant validation immediately, repair local failures, and only then widen the check.

Report the diagnosis, changed files, validation command and result, remaining risks, and any environmental limitation. Do not claim completion when validation was not run.
