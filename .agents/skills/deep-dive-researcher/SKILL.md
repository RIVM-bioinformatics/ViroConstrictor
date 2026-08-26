---
name: deep-dive-researcher
description: Investigate repository behavior, issues, pull requests, reviews, and failures to find root causes and produce evidence-based conclusions. Use when a task asks for a deep investigation or implementation plan.
---

# Deep Dive Researcher

Investigate one concrete question deeply enough to reach a decision. Read `AGENTS.md` and the relevant implementation, tests, history, or review artifacts.

## Procedure

1. Define the question, expected behavior, and the decision the investigation must support.
2. Capture a compact reproduction brief before broad exploration:
	- Trigger: input, command, change, or review comment that exposes the problem.
	- Expected and observed behavior, including the exact error or incorrect result when available.
	- Environment: relevant Python, Snakemake, container, scheduler, platform, and repository versions.
	- Smallest reproducer, or the concrete reason reproduction is unavailable.
	- Evidence already collected and the confidence it supports.
3. Find the nearest code path that computes, mutates, or controls the behavior; do not stop at a wiring layer.
4. Gather the smallest useful evidence set from implementation, callers, tests, configuration, documentation, and version history.
5. Compare the observed behavior with the intended behavior and distinguish confirmed facts from assumptions.
6. Identify edge cases, security implications, compatibility constraints, and affected parallel paths.
7. Choose a cheap discriminating check, such as a focused test, a reproducible command, or a direct call-site comparison.
8. Run the check when execution is requested and use its result to confirm or revise the hypothesis.
9. Produce one conclusion and a concrete next action. Do not repeat an unchanged investigation without new evidence.

## Output

Use this structure when useful:

```markdown
# Investigation

## Conclusion
[Direct answer and confidence]

## Reproduction brief
- Trigger: [input, command, change, or review comment]
- Expected: [intended behavior]
- Observed: [actual behavior or error]
- Environment: [relevant versions and runtime details]
- Reproducer: [smallest reproducer or why unavailable]

## Evidence
- [Implementation or test evidence]
- [Relevant constraint or compatibility fact]

## Risks and edge cases
- [Only material risks]

## Recommended action
[Smallest justified fix, follow-up, or reason no change is needed]

## Validation
[Check run, result, or why execution was unavailable]
```

Keep the brief proportional to the issue. Never fill missing facts with guesses; label assumptions and blockers explicitly.

For review or quality findings, classify each finding as valid, invalid, or requiring clarification, explain whether it conflicts with intended design, and state whether to fix, suppress, or document it. Avoid presenting a plan as completed work.
